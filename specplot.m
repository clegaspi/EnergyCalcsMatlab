%% NQ
% Lets say that we have a ground state energy of theta that looks like
% clear classes
% clf
theta = 0:1:180;
GSBarrier = 10; % kcal/mol
vgs = 0.5*(1-cosd(2*theta));
% figure(1)
% plot(theta,vgs,'r');
%probability as a function of theta
% RT = 0.5922; % I think this may be room temp in kcal/mol
% prob = exp(-vgs/RT);
osc = nq.get_osc();
prob = (nq.esci-nq.esci(1)) .* (osc.^2);
% hold on;
% plot(theta,prob,'b');
% need the excitation energy as a function of theta
% exc = 3.8 + 0.5*(1-cosd(2*theta));
exc = nq.esci-nq.esci(1);
% figure(2)
% plot(theta,exc,'g');
% create a spectrum as a function of energy
ener = 2.5:0.025:7.6; % eV
spec = zeros(size(ener));
LWother = 0.5; % assume a line width of 0.1eV for effects other than torsion
figure(3)
% for it=1:length(theta)
for it=1:length(prob)
   % we assume a gaussian centered at exc(i), with width LWother
   for ie = 1:length(ener)
      g = prob(it) * exp(-(ener-exc(it)).^2./LWother.^2);
      hold on;
      plot(1e7./(ener*8065.73) , g,'r');
      spec = spec + g;
   end
end
% spec = spec./max(spec);
figure(4)
hold on
plot(1e7./(ener*8065.73),spec,'b');

%% DHNQ

% Lets say that we have a ground state energy of theta that looks like
% clear classes
% clf
theta = 0:1:180;
GSBarrier = 10; % kcal/mol
vgs = 0.5*(1-cosd(2*theta));
% figure(1)
% plot(theta,vgs,'r');
%probability as a function of theta
% RT = 0.5922; % I think this may be room temp in kcal/mol
% prob = exp(-vgs/RT);
osc = dhnq.get_osc();
prob = (dhnq.esci-dhnq.esci(1)) .* (osc.^2);
% hold on;
% plot(theta,prob,'b');
% need the excitation energy as a function of theta
% exc = 3.8 + 0.5*(1-cosd(2*theta));
exc = dhnq.esci-dhnq.esci(1);
% figure(2)
% plot(theta,exc,'g');
% create a spectrum as a function of energy
ener = 2.5:0.025:7.6; % eV
spec = zeros(size(ener));
LWother = 0.5; % assume a line width of 0.1eV for effects other than torsion
figure(3)
% for it=1:length(theta)
for it=1:length(prob)
   % we assume a gaussian centered at exc(i), with width LWother
   for ie = 1:length(ener)
      g = prob(it) * exp(-(ener-exc(it)).^2./LWother.^2);
      hold on;
      plot(1e7./(ener*8065.73) , g,'r');
      spec = spec + g;
   end
end
% spec = spec./max(spec);
figure(4)
hold on
plot(1e7./(ener*8065.73),spec,'r');

%% DCNQ

% Lets say that we have a ground state energy of theta that looks like
% clear classes
% clf
theta = 0:1:180;
GSBarrier = 10; % kcal/mol
vgs = 0.5*(1-cosd(2*theta));
% figure(1)
% plot(theta,vgs,'r');
%probability as a function of theta
% RT = 0.5922; % I think this may be room temp in kcal/mol
% prob = exp(-vgs/RT);
osc = dcnq.get_osc();
prob = (dcnq.esci-dcnq.esci(1)) .* (osc.^2);
% hold on;
% plot(theta,prob,'b');
% need the excitation energy as a function of theta
% exc = 3.8 + 0.5*(1-cosd(2*theta));
exc = dcnq.esci-dcnq.esci(1);
% figure(2)
% plot(theta,exc,'g');
% create a spectrum as a function of energy
ener = 2.5:0.025:7.6; % eV
spec = zeros(size(ener));
LWother = 0.5; % assume a line width of 0.1eV for effects other than torsion
figure(3)
% for it=1:length(theta)
for it=1:length(prob)
   % we assume a gaussian centered at exc(i), with width LWother
   for ie = 1:length(ener)
      g = prob(it) * exp(-(ener-exc(it)).^2./LWother.^2);
      hold on;
      plot(1e7./(ener*8065.73) , g,'r');
      spec = spec + g;
   end
end
% spec = spec./max(spec);
figure(4)
plot(1e7./(ener*8065.73),spec,'g');

%% DCDHNQ

% Lets say that we have a ground state energy of theta that looks like
% clear classes
% clf
theta = 0:1:180;
GSBarrier = 10; % kcal/mol
vgs = 0.5*(1-cosd(2*theta));
% figure(1)
% plot(theta,vgs,'r');
%probability as a function of theta
% RT = 0.5922; % I think this may be room temp in kcal/mol
% prob = exp(-vgs/RT);
osc = dcdhnq.get_osc();
prob = (dcdhnq.esci-dcdhnq.esci(1)) .* (osc.^2);
% hold on;
% plot(theta,prob,'b');
% need the excitation energy as a function of theta
% exc = 3.8 + 0.5*(1-cosd(2*theta));
exc = dcdhnq.esci-dcdhnq.esci(1);
% figure(2)
% plot(theta,exc,'g');
% create a spectrum as a function of energy
ener = 2.5:0.025:7.6; % eV
spec = zeros(size(ener));
LWother = 0.5; % assume a line width of 0.1eV for effects other than torsion
figure(3)
% for it=1:length(theta)
for it=1:length(prob)
   % we assume a gaussian centered at exc(i), with width LWother
   for ie = 1:length(ener)
      g = prob(it) * exp(-(ener-exc(it)).^2./LWother.^2);
      hold on;
      plot(1e7./(ener*8065.73) , g,'r');
      spec = spec + g;
   end
end
% spec = spec./max(spec);
figure(4)
plot(1e7./(ener*8065.73),spec,'k');