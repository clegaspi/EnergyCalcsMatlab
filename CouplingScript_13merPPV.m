% %%
% % 13-mer PPV Cross 1
% 
% S = load(['C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\',...
%     'coupling\Exp\NoAngle\13-merPPV-N2Optimized-FrozenGeom-0VA.mat']);
% myexp = S.obj;
% 
% myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?
% 
% efv = myexp.data(1).raw_ampac.r(:,3)-myexp.data(1).raw_ampac.r(:,102);
% 
% mycoup = StateCoupling(2,3,0.0586,0.0592,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_1, fs13_1] = mycoup.run();
% 
% indo13_1 = mycoup.final_indo;
% 
% disp('Crossing 1');
% disp(['Coupling: ', num2str(coup13_1)]);
% disp(['Critical Field: ', num2str(fs13_1)]);
% 
% % % 13-mer PPV Cross 2
% 
% mycoup = StateCoupling(3,4,0.0686,0.0696,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_2, fs13_2] = mycoup.run();
% 
% indo13_2 = mycoup.final_indo;
% 
% disp('Crossing 2');
% disp(['Coupling: ', num2str(coup13_2)]);
% disp(['Critical Field: ', num2str(fs13_2)]);
% 
% % 13-mer PPV Cross 3
% 
% mycoup = StateCoupling(4,5,0.0698,0.0704,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_3, fs13_3] = mycoup.run();
% 
% indo13_3 = mycoup.final_indo;
% 
% disp('Crossing 3');
% disp(['Coupling: ', num2str(coup13_3)]);
% disp(['Critical Field: ', num2str(fs13_3)]);
% 
% % 13-mer PPV Cross 4
% 
% mycoup = StateCoupling(5,6,0.0758,0.0764,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_4, fs13_4] = mycoup.run();
% 
% indo13_4 = mycoup.final_indo;
% 
% disp('Crossing 4');
% disp(['Coupling: ', num2str(coup13_4)]);
% disp(['Critical Field: ', num2str(fs13_4)]);
% 
% % 13-mer PPV Cross 5
% 
% mycoup = StateCoupling(6,7,0.0774,0.0782,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_5, fs13_5] = mycoup.run();
% 
% indo13_5 = mycoup.final_indo;
% 
% disp('Crossing 5');
% disp(['Coupling: ', num2str(coup13_5)]);
% disp(['Critical Field: ', num2str(fs13_5)]);
% 
% % 13-mer PPV Cross 6
% 
% mycoup = StateCoupling(7,8,0.0796,0.0806,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_6, fs13_6] = mycoup.run();
% 
% indo13_6 = mycoup.final_indo;
% 
% disp('Crossing 6');
% disp(['Coupling: ', num2str(coup13_6)]);
% disp(['Critical Field: ', num2str(fs13_6)]);
% 
% % 13-mer PPV Cross 7
% 
% mycoup = StateCoupling(8,9,0.082,0.0828,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_7, fs13_7] = mycoup.run();
% 
% indo13_7 = mycoup.final_indo;
% 
% disp('Crossing 7');
% disp(['Coupling: ', num2str(coup13_7)]);
% disp(['Critical Field: ', num2str(fs13_7)]);
% 
% % 13-mer PPV Cross 8
% 
% mycoup = StateCoupling(9,10,0.0828,0.0832,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_8, fs13_8] = mycoup.run();
% 
% indo13_8 = mycoup.final_indo;
% 
% disp('Crossing 8');
% disp(['Coupling: ', num2str(coup13_8)]);
% disp(['Critical Field: ', num2str(fs13_8)]);
% 
% % 13-mer PPV Cross 9
% 
% mycoup = StateCoupling(10,11,0.0844,0.0848,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_9, fs13_9] = mycoup.run();
% 
% indo13_9 = mycoup.final_indo;
% 
% disp('Crossing 9');
% disp(['Coupling: ', num2str(coup13_9)]);
% disp(['Critical Field: ', num2str(fs13_9)]);
% 
% % 13-mer PPV Cross 10
% 
% mycoup = StateCoupling(11,12,0.0852,0.086,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_10, fs13_10] = mycoup.run();
% 
% indo13_10 = mycoup.final_indo;
% 
% disp('Crossing 10');
% disp(['Coupling: ', num2str(coup13_10)]);
% disp(['Critical Field: ', num2str(fs13_10)]);
% 
% % 13-mer PPV Cross 11
% 
% mycoup = StateCoupling(12,13,0.088,0.0884,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_11, fs13_11] = mycoup.run();
% 
% indo13_11 = mycoup.final_indo;
% 
% disp('Crossing 11');
% disp(['Coupling: ', num2str(coup13_11)]);
% disp(['Critical Field: ', num2str(fs13_11)]);
% 
% % 13-mer PPV Cross 12
% 
% mycoup = StateCoupling(13,14,0.0886,0.09,efv,...
%     'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\13-merPPV\coupling\',...
%     '13-merPPV-N2Optimized');
% 
% [coup13_12, fs13_12] = mycoup.run();
% 
% indo13_12 = mycoup.final_indo;
% 
% disp('Crossing 12');
% disp(['Coupling: ', num2str(coup13_12)]);
% disp(['Critical Field: ', num2str(fs13_12)]);

%% Plots

figure(8);
hold on;

coupling = zeros(1,12);
critfield = zeros(1,12);
states_hdl = zeros(1,12);
ecenter = zeros(1,12);
mingap_energies = zeros(12,2);
mingap_hdl = zeros(1,12);

for i = 1:12
    eval(['coupling(i) = coup13_',num2str(i),';']);
    eval(['critfield(i) = fs13_',num2str(i),';']);
    eval(['tmpindo = indo13_',num2str(i),';']);
    states_hdl(i) = plot(repmat(critfield(i),1,25), tmpindo.esci  - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1), 'm^');
    ecenter(i) = ((tmpindo.esci(2+i) + tmpindo.esci(1+i)) / 2)  - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
    mingap_energies(i,1) = tmpindo.esci(1+i) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
    mingap_energies(i,2) = tmpindo.esci(2+i) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
    mingap_hdl(i) = line(repmat(critfield(i),1,2),mingap_energies(i,:),'Color','r','LineWidth',1.5);
end

coupling_hdl = plot(critfield, ecenter, 'r*');

%% Delete plots

for i = 1:12
    delete(states_hdl(i));
    delete(mingap_hdl(i));
end

delete(coupling_hdl);