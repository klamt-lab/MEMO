load metacyc_universal_model.mat model NADHcons NADPHcons atpm no_carb_met co2ex uptakes acetylcons

%% FVA parameters
cplex_inner= setup_cplex_inner_class_access();
fva_tol= 1e-8;
fvaParams= cplex_inner.ParameterSet.constructor.newInstance([]);
fvaParams.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
fvaParams.setParam(cplex_inner.DoubleParam.EpOpt, fva_tol);
fvaParams.setParam(cplex_inner.DoubleParam.EpRHS, fva_tol);
fvaParams.setParam(cplex_inner.IntParam.RootAlg, 2);

%% select substrates and products
substrates= struct();
substrates.names= {'EX_ACET' 'EX_GLYCEROL' 'EX_ALPHA-GLUCOSE' 'EX_SORBITOL' 'EX_FORMATE'...
  'EX_XYLOSE' 'EX_D-LACTATE' 'EX_METOH' 'EX_SUC'};
[~, substrates.idx]= ismember(substrates.names, model.rxns);
malex= find(strcmp('EX_MAL', model.rxns));

%% additional network configuration
exchange_off= {'EX_HYDROGEN-PEROXIDE'};
[~, exchange_off]= ismember(exchange_off, model.rxns);
reac_off= {'ATPASE-RXN'};
[~, reac_off]= ismember(reac_off, model.rxns);
reac_off= [reac_off NADHcons];

%% thermodynamic parameters
RT= 1.9872036E-3 * (25 + 273.15); % [kc/mol] ;

Cmin= 1e-6*ones(length(model.mets), 1); % mol/l
Cmax= 0.1*ones(length(model.mets), 1); % mol/l
idx= find(strcmp('WATER', model.mets));
Cmin(idx)= 1;
Cmax(idx)= 1;
idx= find(strcmp('PROTON', model.mets));
Cmin(idx)= 1;
Cmax(idx)= 1;
[~, idx]= ismember({'CARBON-DIOXIDE', 'HCO3'}, model.mets);
Cmin(idx)= 1e-7;
Cmax(idx)= 1e-3;

%% cofactor requirement constraints
constr= zeros(0, length(model.rxns));
constr(1, atpm)= 1;
constr(2, NADPHcons)= 1;
constr(3, acetylcons)= 1;
% configuration for 1 ATP, 4 NADPH, 1 acetyl-CoA:
lhs= [1; 4; 1];
rhs= [1; 4; 1];
% alternative configuration for 2 ATP, 3 NADPH, 1 acetyl-CoA:
% lhs= [2; 3; 1];
% rhs= [2; 3; 1];
     
%% structure for storing the results
num_sols= 10;
res= struct();
res(length(substrates.names), 1, num_sols).objval= [];
res(end, end, end).fv= [];
res(end, end, end).status= [];
res(end, end, end).bb= [];
res(end, end, end).fvalb= [];
res(end, end, end).fvaub= [];
res(end, end, end).blocked_fva= [];
res(end, end, end).mmdf= [];

%%
for i= 1:length(substrates.names) % use 'for i= 1' with all substrates open
  i
  j= 1;
  % configure network and run structural FVA
  lplb= model.lb;
  lpub= model.ub;
  lpub(co2ex)= 0; % no net CO2 production
  lpub(malex)= 1000;
  lplb([exchange_off reac_off])= 0;
  lpub([exchange_off reac_off])= 0;
  lplb(setdiff(uptakes, substrates.idx(i)))= 0; % close all uptakes except for the substrates
  % one substrate uptake only:
  lplb(substrates.idx(i))= -10; % limit substrate uptake rate to 10
  % alternative configuration with all substrate uptakes open:
%   lplb(substrates.idx)= -10; % limit substrate uptake rates to 10
  
  fvalhs= [zeros(size(model.S, 1), 1); lhs];
  fvarhs= [zeros(size(model.S, 1), 1); rhs];
  [ind_i, ind_j, val]= find([model.S; constr]);
  fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, fvalhs, fvarhs, lplb, lpub, fvaParams);
  clear ind_i ind_j val;
  if fvaj(3)
    res(i, j, 1).status= 'FVA infeasible';
  else
    res(i, j, 1).fvalb= fvaj(1);
    res(i, j, 1).fvaub= fvaj(2);
    res(i, j, 1).fvalb(abs(res(i, j, 1).fvalb) < fva_tol)= 0;
    res(i, j, 1).fvaub(abs(res(i, j, 1).fvaub) < fva_tol)= 0;

    % setup MILP
    [~, ~, ~, ~, ~, ~, ~, ~, reac_map, obj, z_vars, flux_vars]=...
      max_min_driving_force_pathway(model.S, model.S, [],...
      [], constr, {lhs, rhs}, 1e-7, [], 1000, RT, model.deltaG(:), Cmin, Cmax, [],...
      isnan(model.deltaG(:)), false, res(i, j, 1).fvalb, res(i, j, 1).fvaub, [], 1e-2, NaN, 1);
    to_split= -reac_map(length(model.rxns)+1:end); % backward directions of split reversible reactions
    
    % parameters for the MILP
    obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 80000); % working memory in MB; set as high as your physical RAM allows
    obj.cpx.setParam(cplex_inner.IntParam.FPHeur, 1);
    obj.cpx.setParam(cplex_inner.IntParam.ParallelMode, -1);
    obj.cpx.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.HiddenFeas);
    obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, 3600); % time limit in seconds
    
    % do not count technical reactions in the objective function
    [sel, idx]= ismember([no_carb_met -no_carb_met uptakes -uptakes...
      co2ex -co2ex atpm NADHcons NADPHcons acetylcons], reac_map);
    idx= idx(sel);
    idx= setdiff(1:length(reac_map), idx);
    obj.cpx.getObjective().setExpr(obj.cpx.scalProd(ones(length(reac_map), 1), z_vars(idx)));
    
    for k= 1:num_sols
      obj.cpx.solve(); % run optimization
      
      res(i, j, k).status= char(obj.cpx.getStatus());
      res(i, j, k).bb= obj.cpx.getBestObjValue();
      if obj.cpx.getStatus().equals(cplex_inner.Status.Feasible) || obj.cpx.getStatus().equals(cplex_inner.Status.Optimal)
        res(i, j, k).fv= obj.cpx.getValues(flux_vars);
        res(i, j, k).fv(to_split)= res(i, j, k).fv(to_split) - res(i, j, k).fv(length(model.rxns)+1:end);
        res(i, j, k).fv(length(model.rxns)+1:end)= [];
        res(i, j, k).objval= round(obj.cpx.getObjValue);
        % calculate pathway MDF
        fv= res(i, j, k).fv;
        fv(abs(fv) < 1e-7)= 0;
        res(i, j, k).mmdf= max_min_driving_force_pathway(model.S, model.S, [],...
          [], [], [], 1e-7, [], 1000, RT, model.deltaG, Cmin, Cmax, [],...
          isnan(model.deltaG), false, fv, fv);
        % exclude current solution
        idx= find(fv);
        nr= length(idx);
        [sel, idx]= ismember([idx; -idx], reac_map);
        idx= idx(sel);
        obj.cpx.addLe(obj.cpx.sum(z_vars(idx)), nr - 1);
      elseif obj.cpx.getStatus().equals(cplex_inner.Status.Infeasible)
        break;
      end
    end
  end
  save metacyc_CETCH_screening_multi.mat substrates res
end
