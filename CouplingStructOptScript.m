% path = 'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\coupling\';
% 
% % [xyz, atomtype] = ampac_to_xyz([path,'13-merPPV-Coupling.out']);
% % efv = xyz(3,:)-xyz(102,:);
% % efv = efv / norm(efv);
% 
% 
% S = load([path,'..\13-merPPV\Exp\NoAngle\13-merPPV-0VA.mat']);
% myexp = S.obj;
% myexp.data(1).load_to_memory('indo','load');
% indo = myexp.data(1).raw_indo;
% 
% [rgs, rex, deltar, eszmat, ampac_energy, exindo, ~, numruns] = OptExcStateStructure(path, '13-merPPV-Coupling',...
%     path, '13-merPPV-Coupling',...
%     'field', [0.0 0.0 0.0],...
%     'indo', indo);

