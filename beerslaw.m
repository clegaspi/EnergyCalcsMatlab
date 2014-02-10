figure(2343)
beers = zeros(5,1);
for i=1:5
beers(i) = absdata{i}(993);
end
plot(conc,beers,'bx','MarkerSize',10,'LineWidth',2,'DisplayName','Experimental, \lambda = 304 nm')
m = (conc'\beers);
hold on;
plot(0:1e-5:1.5e-4,m.*[0:1e-5:1.5e-4],'b-','DisplayName', ['Fit, \epsilon \approx ',num2str(round(m)),' M^{-1} cm^{-1}']);
xlabel('Concentration (mol L^{-1})')
ylabel('Absorbance')
title('3,3''-dimethyl-2,2''-Bithiophene - Beer''s Law Fit')
legend('show')

% set(gco, 'DisplayName', 'Experimental, \lambda = 246 nm')
% set(gco, 'DisplayName', 'Fit, \epsilon \approx 7866 M^{-1} cm^{-1}')