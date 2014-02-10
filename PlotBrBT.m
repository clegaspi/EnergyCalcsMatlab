%% Sample B Emission
% Data from 8/20/13
[baxes,babshnd,bem1hnd]=plotyy(babsx,babsy,bem1wl,bem1data);
set(babshnd,'LineWidth',2)
set(bem1hnd,'LineWidth',2,'LineStyle','--')
axes(baxes(2))
hold on
bem2hnd=plot(bem2wl,bem2data,'LineStyle','--','LineWidth',2,'Color','r');
bem3hnd=plot(bem3wl,bem3data,'LineStyle','--','LineWidth',2,'Color','c');
set(baxes(2),'XLimMode','auto')
set(baxes(2),'YLimMode','auto')
set(baxes(2),'YTickMode','auto')
set(baxes(1),'XLim',get(baxes(2),'XLim'))
set(baxes(1),'YLimMode','auto')
set(baxes(1),'YTickMode','auto')
box(baxes(1),'off')
set(bem1hnd,'DisplayName','Emission (\lambda_{ex} = 250 nm)')
set(bem2hnd,'DisplayName','Emission (\lambda_{ex} = 300 nm)')
set(bem3hnd,'DisplayName','Emission (\lambda_{ex} = 400 nm)')
set(babshnd,'DisplayName','Absorbance')
legend('show')
xlabel('Wavelength (nm)')
ylabel('Emission (cps)')
axes(baxes(1))
ylabel('Absorbance')
axes(baxes(2))
title('3,3''-dibromo-2,2''-bithiophene (Sample B) Optical Properties')

%% Sample B Excitation
% Data from 8/20/13
[baxes,babshnd,bex1hnd]=plotyy(babsx,babsy,bex1wl,bex1data);
set(babshnd,'LineWidth',2)
set(bex1hnd,'LineWidth',2,'LineStyle','--')
axes(baxes(2))
hold on
bex2hnd=plot(bex2wl(1:end-30),bex2data(1:end-30),'LineStyle','--','LineWidth',2,'Color','r');
bex3hnd=plot(bex3wl(1:end-40),bex3data(1:end-40),'LineStyle','--','LineWidth',2,'Color','c');
bex4hnd=plot(bex4wl(61:end-20),bex4data(61:end-20),'LineStyle','--','LineWidth',2,'Color','k');
set(baxes(2),'XLimMode','auto')
set(baxes(2),'YLimMode','auto')
set(baxes(2),'YTickMode','auto')
set(baxes(1),'XLim',get(baxes(2),'XLim'))
set(baxes(1),'YLimMode','auto')
set(baxes(1),'YTickMode','auto')
box(baxes(1),'off')
set(bex1hnd,'DisplayName','Excitation (\lambda_{em} = 290 nm)')
set(bex2hnd,'DisplayName','Excitation (\lambda_{em} = 390 nm)')
set(bex3hnd,'DisplayName','Excitation (\lambda_{em} = 430 nm)')
set(bex4hnd,'DisplayName','Excitation (\lambda_{em} = 500 nm)')
set(babshnd,'DisplayName','Absorbance')
legend('show')
xlabel('Wavelength (nm)')
ylabel('Emission (cps)')
axes(baxes(1))
ylabel('Absorbance')
axes(baxes(2))
title('3,3''-dibromo-2,2''-bithiophene (Sample B) Optical Properties')

%% Sample D Emission
% Data from 8/23/13
[daxes,dabshnd,dem1hnd]=plotyy(dabsx,dabsy,dem1wl,dem1data);
set(dabshnd,'LineWidth',2)
set(dem1hnd,'LineWidth',2,'LineStyle','--')
axes(daxes(2))
hold on
dem2hnd=plot(dem2wl,dem2data,'LineStyle','--','LineWidth',2,'Color','r');
dem3hnd=plot(dem3wl,dem3data,'LineStyle','--','LineWidth',2,'Color','c');
dem4hnd=plot(dem4wl,dem4data,'LineWidth',2,'LineStyle','--','Color','k');
set(daxes(2),'XLimMode','auto')
set(daxes(2),'YLimMode','auto')
set(daxes(2),'YTickMode','auto')
set(daxes(1),'XLim',get(daxes(2),'XLim'))
set(daxes(1),'YLimMode','auto')
set(daxes(1),'YTickMode','auto')
box(daxes(1),'off')
set(dem1hnd,'DisplayName','Emission (\lambda_{ex} = 240 nm)')
set(dem2hnd,'DisplayName','Emission (\lambda_{ex} = 280 nm)')
set(dem3hnd,'DisplayName','Emission (\lambda_{ex} = 320 nm)')
set(dem4hnd,'DisplayName','Emission (\lambda_{ex} = 400 nm)')
set(dabshnd,'DisplayName','Absorbance')
legend('show')
xlabel('Wavelength (nm)')
ylabel('Emission (cps)')
axes(daxes(1))
ylabel('Absorbance')
axes(daxes(2))
title('5,5''-dibromo-2,2''-bithiophene (Sample D) Optical Properties')

%% Sample D Excitation
% Data from 8/23/13
[daxes,dabshnd,dex1hnd]=plotyy(dabsx,dabsy,dex1wl,dex1data);
set(dabshnd,'LineWidth',2)
set(dex1hnd,'LineWidth',2,'LineStyle','--')
axes(daxes(2))
hold on
dex2hnd=plot(dex2wl,dex2data,'LineStyle','--','LineWidth',2,'Color','r');
dex3hnd=plot(dex3wl,dex3data,'LineStyle','--','LineWidth',2,'Color','c');
dex4hnd=plot(dex4wl(60:end),dex4data(60:end),'LineWidth',2,'LineStyle','--','Color','k');
set(daxes(2),'XLimMode','auto')
set(daxes(2),'YLimMode','auto')
set(daxes(2),'YTickMode','auto')
set(daxes(1),'XLim',get(daxes(2),'XLim'))
set(daxes(1),'YLimMode','auto')
set(daxes(1),'YTickMode','auto')
box(daxes(1),'off')
set(dex1hnd,'DisplayName','Excitation (\lambda_{em} = 320 nm)')
set(dex2hnd,'DisplayName','Excitation (\lambda_{em} = 370 nm)')
set(dex3hnd,'DisplayName','Excitation (\lambda_{em} = 390 nm)')
set(dex4hnd,'DisplayName','Excitation (\lambda_{em} = 500 nm)')
set(dabshnd,'DisplayName','Absorbance')
legend('show')
xlabel('Wavelength (nm)')
ylabel('Emission (cps)')
axes(daxes(1))
ylabel('Absorbance')
axes(daxes(2))
title('5,5''-dibromo-2,2''-bithiophene (Sample D) Optical Properties')