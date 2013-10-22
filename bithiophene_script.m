%% Bithiophene
twist = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\trans-bithiophene\trans-bithiophene-slighttwist.ido');
planar = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\trans-bithiophene\trans-bithiophene-planar.ido');

%% Terthiophene
twist = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-terthiophene\alltrans-terthiophene-twist.ido');
planar = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-terthiophene\alltrans-terthiophene-planar.ido');

%% Quaterthiophene
twist = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-quaterthiophene\alltrans-quaterthiophene-twist.ido');
planar = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-quaterthiophene\alltrans-quaterthiophene-planar.ido');

%% Quinquithiophene
twist = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-quinquithiophene\alltrans-quinquithiophene-twist.ido');
planar = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-quinquithiophene\alltrans-quinquithiophene-planar.ido');

%% Sexithiophene
twist = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-sexithiophene\alltrans-sexithiophene-twist.ido');
planar = Indo.LoadExistingData('C:\Users\clegaspi\Documents\MATLAB\data\alltrans-sexithiophene\alltrans-sexithiophene-planar.ido');

%% Plot and calculate

fignum = 2;
figure(fignum)
hold off
plot(zeros(1,25),twist.esci, 'g^')
hold on
twist_os = twist.get_osc();
twist_os = twist_os .* 30 ./ max(twist_os) + 1e-3;
for i = 1:25
    plot(0,twist.esci(i), 'bo', 'MarkerSize',twist_os(i))
end
figure(fignum);
plot(ones(1,25),planar.esci, 'g^')
planar_os = planar.get_osc();
planar_os = planar_os .* 30 ./ max(planar_os) + 1e-3;
for i = 1:25
    plot(1,planar.esci(i), 'bo', 'MarkerSize',planar_os(i))
end

stabilization = (twist.esci(2)-twist.esci(1)) - (planar.esci(2)-planar.esci(1));
disp(['Stabilization: ',num2str(stabilization)]);
