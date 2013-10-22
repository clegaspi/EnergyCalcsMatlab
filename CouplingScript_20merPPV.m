% S = load(['C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\',...
%     'coupling\Exp\NoAngle\20-merPPV-N2Optimized-0VA.mat']);
% myexp = S.obj;
% 
% myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?
% 
% efv = myexp.data(1).raw_ampac.r(:,3)-myexp.data(1).raw_ampac.r(:,158);
% 
% mycoup = StateCoupling(2,3,0.0330,0.0338,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
%     '20-merPPV-N2Optimized');
% 
% [coup20_1, fs20_1] = mycoup.run();
% 
% indo20_1 = mycoup.final_indo;
% 
% disp('Crossing 1');
% disp(['Coupling: ', num2str(coup20_1)]);
% disp(['Critical Field: ', num2str(fs20_1)]);
% 
% mycoup = StateCoupling(3,4,0.0374,0.0380,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
%     '20-merPPV-N2Optimized');
% 
% [coup20_2, fs20_2] = mycoup.run();
% 
% indo20_2 = mycoup.final_indo;
% 
% disp('Crossing 2');
% disp(['Coupling: ', num2str(coup20_2)]);
% disp(['Critical Field: ', num2str(fs20_2)]);
% 
% mycoup = StateCoupling(4,5,0.0380,0.0384,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
%     '20-merPPV-N2Optimized');
% 
% [coup20_3, fs20_3] = mycoup.run();
% 
% indo20_3 = mycoup.final_indo;
% 
% disp('Crossing 3');
% disp(['Coupling: ', num2str(coup20_3)]);
% disp(['Critical Field: ', num2str(fs20_3)]);
% 
% mycoup = StateCoupling(5,6,0.0412,0.0418,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
%     '20-merPPV-N2Optimized');
% 
% [coup20_4, fs20_4] = mycoup.run();
% 
% indo20_4 = mycoup.final_indo;
% 
% disp('Crossing 4');
% disp(['Coupling: ', num2str(coup20_4)]);
% disp(['Critical Field: ', num2str(fs20_4)]);
% 
% mycoup = StateCoupling(6,7,0.0422,0.0426,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
%     '20-merPPV-N2Optimized');
% 
% [coup20_5, fs20_5] = mycoup.run();
% 
% indo20_5 = mycoup.final_indo;
% 
% disp('Crossing 5');
% disp(['Coupling: ', num2str(coup20_5)]);
% disp(['Critical Field: ', num2str(fs20_5)]);
%%
mycoup = StateCoupling(7,8,0.0430,0.0434,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
    '20-merPPV-N2Optimized');

[coup20_6, fs20_6] = mycoup.run();

indo20_6 = mycoup.final_indo;

disp('Crossing 6');
disp(['Coupling: ', num2str(coup20_6)]);
disp(['Critical Field: ', num2str(fs20_6)]);

mycoup = StateCoupling(8,9,0.0449,0.0450,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
    '20-merPPV-N2Optimized');

[coup20_7, fs20_7] = mycoup.run();

indo20_7 = mycoup.final_indo;

disp('Crossing 7');
disp(['Coupling: ', num2str(coup20_7)]);
disp(['Critical Field: ', num2str(fs20_7)]);

mycoup = StateCoupling(9,10,0.0450,0.0451,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
    '20-merPPV-N2Optimized');

[coup20_8, fs20_8] = mycoup.run();

indo20_8 = mycoup.final_indo;

disp('Crossing 8');
disp(['Coupling: ', num2str(coup20_8)]);
disp(['Critical Field: ', num2str(fs20_8)]);

mycoup = StateCoupling(10,11,0.0462,0.0466,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\20-merPPV\coupling\',...
    '20-merPPV-N2Optimized');

[coup20_9, fs20_9] = mycoup.run();

indo20_9 = mycoup.final_indo;

disp('Crossing 9');
disp(['Coupling: ', num2str(coup20_9)]);
disp(['Critical Field: ', num2str(fs20_9)]);
%% Plots

ncross = 9;
nstates = 12;

figure(8);
hold on;

coupling = zeros(1,ncross);
critfield = zeros(1,ncross);
states_hdl = zeros(1,ncross);
ecenter = zeros(1,ncross);
mingap_energies = zeros(ncross,2);
mingap_hdl = zeros(1,ncross);

for i = 1:ncross
    eval(['coupling(i) = coup20_',num2str(i),';']);
    eval(['critfield(i) = fs20_',num2str(i),';']);
    eval(['tmpindo = indo20_',num2str(i),';']);
    states_hdl(i) = plot(repmat(critfield(i),1,nstates), tmpindo.esci, 'm^');
    
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