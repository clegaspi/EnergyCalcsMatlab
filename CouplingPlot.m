%% Load data

% load('c:\Users\clegaspi\Documents\MATLAB\data\coupling-PPV-relaxandnon.mat');
load('c:\Users\clegaspi\Documents\MATLAB\data\coupling-alloligomers-relaxandnon.mat');

%% log Coupling vs. Crit Field
% plot(crit_8_relax,log10(coup_8_relax),'b*')
% hold on
% plot(crit_13_relax_nooutliers,log10(coup_13_relax_nooutliers),'ro')
% plot(crit_15_relax_nooutliers,log10(coup_15_relax_nooutliers),'g+')
% plot(crit_17_relax_nooutliers,log10(coup_17_relax_nooutliers),'mx')
% plot(crit_20_relax_nooutliers,log10(coup_20_relax_nooutliers),'kd')
% plot(crit_6to20_norelax, log10(coup_6to20_norelax), 'c^')

plot(crit_8_relax,log10(coup_8_relax),'b-','DisplayName','OPPV-8 (Relaxed)')
hold on
plot(crit_13_relax_nooutliers,log10(coup_13_relax_nooutliers),'r-','DisplayName','OPPV-13 (Relaxed)')
plot(crit_15_relax_nooutliers,log10(coup_15_relax_nooutliers),'g-','DisplayName','OPPV-15 (Relaxed)')
plot(crit_17_relax_nooutliers,log10(coup_17_relax_nooutliers),'m-','DisplayName','OPPV-17 (Relaxed)')
plot(crit_20_relax_nooutliers,log10(coup_20_relax_nooutliers),'k-','DisplayName','OPPV-20 (Relaxed)')
plot(crit_6to20_norelax, log10(coup_6to20_norelax), 'c-','DisplayName',...
    sprintf('%s\n%s','OPPV-6 to 20 (Unrelaxed)','First Crossing Only'))

plot(pfh8_fs_nooutliers,log10(pfh8_coup_nooutliers),'bo','DisplayName','PF-8 (Relaxed)')
plot(pfh13_fs_nooutliers,log10(pfh13_coup_nooutliers),'ro','DisplayName','PF-13 (Relaxed)')
plot(melppp8_fs_nooutliers,log10(melppp8_coup_nooutliers),'b*','DisplayName','LP-8 (Relaxed)')
plot(melppp13_fs_nooutliers,log10(melppp13_coup_nooutliers),'r*','DisplayName','LP-13 (Relaxed)')

xlabel('Critical Field / [V cm^-^1] x 10^-^8')
ylabel('log10( Coupling Energy ) / log10( [eV] )')
title('Coupling Energy of Crossing States vs. Critical Field of Crossing for Various Conjugated Oligomers')
legend show

%% Coupling vs. Crit Field
% plot(crit_8_relax,coup_8_relax,'b*')
% hold on
% plot(crit_13_relax_nooutliers,coup_13_relax_nooutliers,'ro')
% plot(crit_15_relax_nooutliers,coup_15_relax_nooutliers,'g+')
% plot(crit_17_relax_nooutliers,coup_17_relax_nooutliers,'mx')
% plot(crit_20_relax_nooutliers,coup_20_relax_nooutliers,'kd')
% plot(crit_6to20_norelax, coup_6to20_norelax, 'c^')
% 
% xlabel('Critical Field / [V cm^-^1] x 10^-^8')
% ylabel('Coupling Energy / [eV]')
% title('Coupling Energy of Crossing States vs. Critical Field of Crossing for Various OPPV Oligomers')
% legend('OPPV-8 (Relaxed)','OPPV-13 (Relaxed)','OPPV-15 (Relaxed)','OPPV-17 (Relaxed)','OPPV-20 (Relaxed)',...
%     sprintf('%s\n%s','OPPV-6 to 20 (Unrelaxed)','First Crossing Only'))

plot(crit_8_relax,coup_8_relax,'b-','DisplayName','OPPV-8 (Relaxed)')
hold on
plot(crit_13_relax_nooutliers,coup_13_relax_nooutliers,'r-','DisplayName','OPPV-13 (Relaxed)')
plot(crit_15_relax_nooutliers,coup_15_relax_nooutliers,'g-','DisplayName','OPPV-15 (Relaxed)')
plot(crit_17_relax_nooutliers,coup_17_relax_nooutliers,'m-','DisplayName','OPPV-17 (Relaxed)')
plot(crit_20_relax_nooutliers,coup_20_relax_nooutliers,'k-','DisplayName','OPPV-20 (Relaxed)')
plot(crit_6to20_norelax, coup_6to20_norelax, 'c-','DisplayName',...
    sprintf('%s\n%s','OPPV-6 to 20 (Unrelaxed)','First Crossing Only'))

plot(pfh8_fs_nooutliers,pfh8_coup_nooutliers,'bo','DisplayName','PF-8 (Relaxed)')
plot(pfh13_fs_nooutliers,pfh13_coup_nooutliers,'ro','DisplayName','PF-13 (Relaxed)')
plot(melppp8_fs_nooutliers,melppp8_coup_nooutliers,'b*','DisplayName','LP-8 (Relaxed)')
plot(melppp13_fs_nooutliers,melppp13_coup_nooutliers,'r*','DisplayName','LP-13 (Relaxed)')

xlabel('Critical Field / [V cm^-^1] x 10^-^8')
ylabel('Coupling Energy / [eV]')
title('Coupling Energy of Crossing States vs. Critical Field of Crossing for Various Conjugated Oligomers')
legend show

