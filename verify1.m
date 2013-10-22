% %% test of coulomb function
% clear classes;
% t1 = Tchain(50);
% t1.U = 0.0; t1.eps = 2.0; t1.betaE = -1.0;
% r = 0:0.1:20;
% coul = t1.hubbard(r);
% plot(r,coul,'b');
% hold on;
% plot(r(10:end),14.397./(t1.eps * r(10:end)),'r')
% 
% %% what is the bandwidth
% a = zeros(100,100);
% for i=1:99
%    a(i,i+1) = -1;
%    a(i+1,i) = -1;
% end
% [ev,ener] = eig(a);
% %%
% close all;
% clear classes;
% t1 = Tchain(21,1); % use 2n+1 sites
% t1.U = 0.0; t1.eps = 2.0; t1.betaH = 1.0; t1.betaE = -1.0;
% t1.aSite = 7; t1.Eo = 7.5;
% field = 0.0;
% 
% ic=0;
% for U = 0:0.5:2
%    t1.U = U;
%    ham = t1.sciMatrix(field);
%    [ev,estupid] = eig(ham);
%    ener = diag(estupid);
%    tm = ev' * (t1.transVector)';
%    maxTM = max(abs(tm));
%    ic = ic+1;
%    sumrule(ic) = sum(ener .* tm.^2);
%    figure(1);
%    hold on;
%    for i=1:t1.nbasis
%       plot(U,ener(i),'b.');
%       plot(U,ener(i),'ro','markerSize',(abs(tm(i))*7/maxTM) +1e-3);
%    end
%    prop = ev(:,1).^2;
%    ehsep = t1.EHDiff;
%    psep = zeros(t1.nsites+1,1);
%    for isep = 0:t1.nsites
%       psep(isep+1) = sum(prop(ehsep == isep));
%    end
%    psep(1) = 2.0 * psep(1);
%    figure(2)
%    hold on;
%    plot(0:t1.nsites,psep,'r-o');
% end
%%
close all;
% clear classes;
t1 = Tchain(7,0); % use 2n+1 sites
t1.U = 1.5; t1.eps = 2.0; t1.betaH = 1.0; t1.betaE = -1.0;
t1.aSite = 7; t1.Eo = 7.5;
field = 0.0;

ic=0;
for field = 0:0.001:0.02
   ham = t1.sciMatrix(field);
   [ev,estupid] = eig(ham);
   ener = diag(estupid);
   tm = ev' * (t1.transVector)';
   maxTM = max(abs(tm));
   ic = ic+1;
   sumrule(ic) = sum(ener .* tm.^2);
   figure(1);
   hold on;
   for i=1:t1.nbasis
      plot(field,ener(i),'b.');
      plot(field,ener(i),'ro','markerSize',(abs(tm(i))*7/maxTM) +1e-3);
   end
   prop = ev(:,1).^2;
   ehsep = t1.EHDiff;
   psep = zeros(t1.nsites+1,1);
   for isep = 0:t1.nsites
      psep(isep+1) = sum(prop(ehsep == isep));
   end
   psep(1) = 2.0 * psep(1);
   figure(2)
   hold on;
   plot(0:t1.nsites,psep,'r-o');
end




