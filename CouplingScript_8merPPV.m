S = load(['C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\8-merPPV\',...
    'coupling\Exp\NoAngle\8-merPPV-N2Optimized-0VA.mat']);
myexp = S.obj;

myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?

efv = myexp.data(1).raw_ampac.r(:,6)-myexp.data(1).raw_ampac.r(:,62);

mycoup = StateCoupling(2,3,0.11,0.12,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\8-merPPV\coupling\',...
    '8-merPPV-N2Optimized');

[coup8_1, fs8_1] = mycoup.run();

indo8_1 = mycoup.final_indo;

disp('Crossing 1');
disp(['Coupling: ', num2str(coup8_1)]);
disp(['Critical Field: ', num2str(fs8_1)]);

mycoup = StateCoupling(3,6,0.14,0.15,efv,...
    'C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\8-merPPV\coupling\',...
    '8-merPPV-N2Optimized');

[coup8_2, fs8_2] = mycoup.run();

indo8_2 = mycoup.final_indo;

disp('Crossing 2');
disp(['Coupling: ', num2str(coup8_2)]);
disp(['Critical Field: ', num2str(fs8_2)]);

%% Plots

figure(8);
hold on;

coupling = zeros(1,2);
critfield = zeros(1,2);
states_hdl = zeros(1,2);
ecenter = zeros(1,2);
mingap_energies = zeros(2,2);
mingap_hdl = zeros(1,2);

for i = 1:2
    eval(['coupling(i) = coup8_',num2str(i),';']);
    eval(['critfield(i) = fs8_',num2str(i),';']);
    eval(['tmpindo = indo8_',num2str(i),';']);
    states_hdl(i) = plot(repmat(critfield(i),1,25), tmpindo.esci - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1), 'm^');
    if (i == 1)
        ecenter(i) = ((tmpindo.esci(2) + tmpindo.esci(3)) / 2) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
        mingap_energies(i,1) = tmpindo.esci(2) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
        mingap_energies(i,2) = tmpindo.esci(3) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
    else
        ecenter(i) = ((tmpindo.esci(3) + tmpindo.esci(6)) / 2) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
        mingap_energies(i,1) = tmpindo.esci(3) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
        mingap_energies(i,2) = tmpindo.esci(6) - tmpindo.esci(1) + tmpindo.hfE - igs{1,1}(1);
    end
    mingap_hdl(i) = line(repmat(critfield(i),1,2),mingap_energies(i,:),'Color','r','LineWidth',1.5);
end

coupling_hdl = plot(critfield, ecenter, 'r*');

%% Delete plots

for i = 1:12
    delete(states_hdl(i));
    delete(mingap_hdl(i));
end

delete(coupling_hdl);