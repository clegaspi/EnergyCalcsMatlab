% S = load(['C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\',...
%     'coupling\Exp\NoAngle\15-merPPV-N2Optimized-0VA.mat']);
% myexp = S.obj;
% 
% myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?
% 
% efv = myexp.data(1).raw_ampac.r(:,3)-myexp.data(1).raw_ampac.r(:,118);
% 
% mycoup = StateCoupling(2,3,0.047,0.048,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_1, fs15_1] = mycoup.run();
% 
% indo15_1 = mycoup.final_indo;
% 
% disp('Crossing 1');
% disp(['Coupling: ', num2str(coup15_1)]);
% disp(['Critical Field: ', num2str(fs15_1)]);
% 
% mycoup = StateCoupling(3,4,0.055,0.0552,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_2, fs15_2] = mycoup.run();
% 
% indo15_2 = mycoup.final_indo;
% 
% disp('Crossing 2');
% disp(['Coupling: ', num2str(coup15_2)]);
% disp(['Critical Field: ', num2str(fs15_2)]);
% 
% mycoup = StateCoupling(4,5,0.0558,0.056,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_3, fs15_3] = mycoup.run();
% 
% indo15_3 = mycoup.final_indo;
% 
% disp('Crossing 3');
% disp(['Coupling: ', num2str(coup15_3)]);
% disp(['Critical Field: ', num2str(fs15_3)]);
% 
% mycoup = StateCoupling(5,6,0.0616,0.0618,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_4, fs15_4] = mycoup.run();
% 
% indo15_4 = mycoup.final_indo;
% 
% disp('Crossing 4');
% disp(['Coupling: ', num2str(coup15_4)]);
% disp(['Critical Field: ', num2str(fs15_4)]);
% 
% mycoup = StateCoupling(6,7,0.0628,0.063,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_5, fs15_5] = mycoup.run();
% 
% indo15_5 = mycoup.final_indo;
% 
% disp('Crossing 5');
% disp(['Coupling: ', num2str(coup15_5)]);
% disp(['Critical Field: ', num2str(fs15_5)]);
% 
% mycoup = StateCoupling(7,8,0.0634,0.0636,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_6, fs15_6] = mycoup.run();
% 
% indo15_6 = mycoup.final_indo;
% 
% disp('Crossing 6');
% disp(['Coupling: ', num2str(coup15_6)]);
% disp(['Critical Field: ', num2str(fs15_6)]);
% 
% mycoup = StateCoupling(8,9,0.0648,0.0650,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_7, fs15_7] = mycoup.run();
% 
% indo15_7 = mycoup.final_indo;
% 
% disp('Crossing 7');
% disp(['Coupling: ', num2str(coup15_7)]);
% disp(['Critical Field: ', num2str(fs15_7)]);
% 
% mycoup = StateCoupling(9,10,0.0672,0.0676,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_8, fs15_8] = mycoup.run();
% 
% indo15_8 = mycoup.final_indo;
% 
% disp('Crossing 8');
% disp(['Coupling: ', num2str(coup15_8)]);
% disp(['Critical Field: ', num2str(fs15_8)]);
% 
% mycoup = StateCoupling(10,11,0.0676,0.0678,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_9, fs15_9] = mycoup.run();
% 
% indo15_9 = mycoup.final_indo;
% 
% disp('Crossing 9');
% disp(['Coupling: ', num2str(coup15_9)]);
% disp(['Critical Field: ', num2str(fs15_9)]);
% 
% mycoup = StateCoupling(11,12,0.0692,0.0694,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_10, fs15_10] = mycoup.run();
% 
% indo15_10 = mycoup.final_indo;
% 
% disp('Crossing 10');
% disp(['Coupling: ', num2str(coup15_10)]);
% disp(['Critical Field: ', num2str(fs15_10)]);

mycoup = StateCoupling(12,13,0.0702,0.0706,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
    '15-merPPV-N2Optimized');

[coup15_11, fs15_11] = mycoup.run();

indo15_11 = mycoup.final_indo;

disp('Crossing 11');
disp(['Coupling: ', num2str(coup15_11)]);
disp(['Critical Field: ', num2str(fs15_11)]);

mycoup = StateCoupling(13,14,0.0718,0.072,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
    '15-merPPV-N2Optimized');

[coup15_12, fs15_12] = mycoup.run();

indo15_12 = mycoup.final_indo;

disp('Crossing 12');
disp(['Coupling: ', num2str(coup15_12)]);
disp(['Critical Field: ', num2str(fs15_12)]);

% mycoup = StateCoupling(14,15,?,?,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\15-merPPV\coupling\',...
%     '15-merPPV-N2Optimized');
% 
% [coup15_13, fs15_13] = mycoup.run();
% 
% indo15_13 = mycoup.final_indo;
% 
% disp('Crossing 13');
% disp(['Coupling: ', num2str(coup15_13)]);
% disp(['Critical Field: ', num2str(fs15_13)]);
%% Plots

ncross = 12;

figure(8);
hold on;

coupling = zeros(1,ncross);
critfield = zeros(1,ncross);
states_hdl = zeros(1,ncross);
ecenter = zeros(1,ncross);
mingap_energies = zeros(ncross,2);
mingap_hdl = zeros(1,ncross);

for i = 1:ncross
    eval(['coupling(i) = coup15_',num2str(i),';']);
    eval(['critfield(i) = fs15_',num2str(i),';']);
    eval(['tmpindo = indo15_',num2str(i),';']);
    states_hdl(i) = plot(repmat(critfield(i),1,length(tmpindo.esci)), tmpindo.esci, 'm^');
    
    ecenter(i) = (tmpindo.esci(i+2) + tmpindo.esci(i+1)) / 2;
    mingap_energies(i,1) = tmpindo.esci(i+1);
    mingap_energies(i,2) = tmpindo.esci(i+2);
    mingap_hdl(i) = line(repmat(critfield(i),1,2),mingap_energies(i,:),'Color','r','LineWidth',1.5);
end

coupling_hdl = plot(critfield, ecenter, 'r*');

%% Delete plots

for i = 1:ncross
    delete(states_hdl(i));
    delete(mingap_hdl(i));
end

delete(coupling_hdl);