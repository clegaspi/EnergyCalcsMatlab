% %% 3-mer PPV
% 
% [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\3-merPPV-Coupling.out');
% efv = xyz(3,:)-xyz(22,:);
% 
% mycoup = StateCoupling(0.36,1.0,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '3-merPPV-Coupling');
% 
% [coup3, fs3] = mycoup.run();
% % 
% % %% 6-mer PPV
% % 
% % [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\6-merPPV-Coupling.out');
% % efv = xyz(3,:)-xyz(46,:);
% % 
% % mycoup = StateCoupling(0.145,0.157,efv,...
% %     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
% %     '6-merPPV-Coupling');
% % 
% % [coup6, fs6] = mycoup.run();
% % 
% % %% 7-mer PPV
% % 
% % [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\7-merPPV-Coupling.out');
% % efv = xyz(3,:)-xyz(54,:);
% % 
% % mycoup = StateCoupling(0.115,0.127,efv,...
% %     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
% %     '7-merPPV-Coupling');
% % 
% % [coup7, fs7] = mycoup.run();
% 
% %% 8-mer PPV
% [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\8-merPPV-Coupling.out');
% efv = xyz(6,:)-xyz(62,:);
% 
% mycoup = StateCoupling(0.095,0.105,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '8-merPPV-Coupling');
% 
% [coup8, fs8] = mycoup.run();
% 
% 
% %% 9-mer PPV
% 
% [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\9-merPPV-Coupling.out');
% efv = xyz(3,:)-xyz(70,:);
% 
% mycoup = StateCoupling(0.079,0.091,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '9-merPPV-Coupling');
% 
% [coup9, fs9] = mycoup.run();
% 
% %% 10-mer PPV
% 
% [xyz, atomtype] = ampac_to_xyz('C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\10-merPPV-Coupling.out');
% efv = xyz(3,:)-xyz(78,:);
% 
% mycoup = StateCoupling(0.069,0.080,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '10-merPPV-Coupling');
% 
% [coup10, fs10] = mycoup.run();

