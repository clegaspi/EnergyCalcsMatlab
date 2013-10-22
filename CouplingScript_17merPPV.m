% S = load(['C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\',...
%     'coupling\Exp\NoAngle\17-merPPV-N2Optimized-0VA.mat']);
% myexp = S.obj;
% 
% myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?
% 
% efv = myexp.data(1).raw_ampac.r(:,3)-myexp.data(1).raw_ampac.r(:,134);
% 
% mycoup = StateCoupling(2,3,0.04,0.041,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_1, fs17_1] = mycoup.run();
% 
% indo17_1 = mycoup.final_indo;
% 
% disp('Crossing 1');
% disp(['Coupling: ', num2str(coup17_1)]);
% disp(['Critical Field: ', num2str(fs17_1)]);
% 
% mycoup = StateCoupling(3,4,0.046,0.047,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_2, fs17_2] = mycoup.run();
% 
% indo17_2 = mycoup.final_indo;
% 
% disp('Crossing 2');
% disp(['Coupling: ', num2str(coup17_2)]);
% disp(['Critical Field: ', num2str(fs17_2)]);
% 
% mycoup = StateCoupling(4,5,0.047,0.048,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_3, fs17_3] = mycoup.run();
% 
% indo17_3 = mycoup.final_indo;
% 
% disp('Crossing 3');
% disp(['Coupling: ', num2str(coup17_3)]);
% disp(['Critical Field: ', num2str(fs17_3)]);
% 
% mycoup = StateCoupling(5,6,0.0516,0.0520,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_4, fs17_4] = mycoup.run();
% 
% indo17_4 = mycoup.final_indo;
% 
% disp('Crossing 4');
% disp(['Coupling: ', num2str(coup17_4)]);
% disp(['Critical Field: ', num2str(fs17_4)]);
% 
% mycoup = StateCoupling(6,7,0.053,0.0532,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_5, fs17_5] = mycoup.run();
% 
% indo17_5 = mycoup.final_indo;
% 
% disp('Crossing 5');
% disp(['Coupling: ', num2str(coup17_5)]);
% disp(['Critical Field: ', num2str(fs17_5)]);
% 
% mycoup = StateCoupling(7,8,0.0540,0.05422,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_6, fs17_6] = mycoup.run();
% 
% indo17_6 = mycoup.final_indo;
% 
% disp('Crossing 6');
% disp(['Coupling: ', num2str(coup17_6)]);
% disp(['Critical Field: ', num2str(fs17_6)]);
% 
% mycoup = StateCoupling(8,9,0.05422,0.0544,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_7, fs17_7] = mycoup.run();
% 
% indo17_7 = mycoup.final_indo;
% 
% disp('Crossing 7');
% disp(['Coupling: ', num2str(coup17_7)]);
% disp(['Critical Field: ', num2str(fs17_7)]);
% 
% mycoup = StateCoupling(9,10,0.0562,0.0566,efv,...
%     'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
%     '17-merPPV-N2Optimized');
% 
% [coup17_8, fs17_8] = mycoup.run();
% 
% indo17_8 = mycoup.final_indo;
% 
% disp('Crossing 8');
% disp(['Coupling: ', num2str(coup17_8)]);
% disp(['Critical Field: ', num2str(fs17_8)]);

mycoup = StateCoupling(10,11,0.0576,0.058,efv,...
    'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
    '17-merPPV-N2Optimized');

[coup17_9, fs17_9] = mycoup.run();

indo17_9 = mycoup.final_indo;

disp('Crossing 9');
disp(['Coupling: ', num2str(coup17_9)]);
disp(['Critical Field: ', num2str(fs17_9)]);

mycoup = StateCoupling(11,12,0.0584,0.0588,efv,...
    'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
    '17-merPPV-N2Optimized');

[coup17_10, fs17_10] = mycoup.run();

indo17_10 = mycoup.final_indo;

disp('Crossing 10');
disp(['Coupling: ', num2str(coup17_10)]);
disp(['Critical Field: ', num2str(fs17_10)]);

mycoup = StateCoupling(12,13,0.0598,0.06,efv,...
    'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
    '17-merPPV-N2Optimized');

[coup17_11, fs17_11] = mycoup.run();

indo17_11 = mycoup.final_indo;

disp('Crossing 11');
disp(['Coupling: ', num2str(coup17_11)]);
disp(['Critical Field: ', num2str(fs17_11)]);

mycoup = StateCoupling(13,14,0.06,0.0604,efv,...
    'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
    '17-merPPV-N2Optimized');

[coup17_12, fs17_12] = mycoup.run();

indo17_12 = mycoup.final_indo;

disp('Crossing 12');
disp(['Coupling: ', num2str(coup17_12)]);
disp(['Critical Field: ', num2str(fs17_12)]);

mycoup = StateCoupling(14,15,0.0604,0.0608,efv,...
    'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\coupling\',...
    '17-merPPV-N2Optimized');

[coup17_13, fs17_13] = mycoup.run();

indo17_13 = mycoup.final_indo;

disp('Crossing 13');
disp(['Coupling: ', num2str(coup17_13)]);
disp(['Critical Field: ', num2str(fs17_13)]);
%% Plots

ncross = 13;

figure(8);
hold on;

coupling = zeros(1,ncross);
critfield = zeros(1,ncross);
states_hdl = zeros(1,ncross);
ecenter = zeros(1,ncross);
mingap_energies = zeros(ncross,2);
mingap_hdl = zeros(1,ncross);

for i = 1:ncross
    eval(['coupling(i) = coup17_',num2str(i),';']);
    eval(['critfield(i) = fs17_',num2str(i),';']);
    eval(['tmpindo = indo17_',num2str(i),';']);
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