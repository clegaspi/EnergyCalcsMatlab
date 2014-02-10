%%
emwlval = [254 304 347];
exwlval = [374 424];
plottitle = '2TM - 1x';

%%
emhdl = zeros(length(emdata),1);
color = 'rgkmc';
figure(52)
[axhdl,abshdl,emhdl(1)]=plotyy(abswl{1},absdata{1},emwl{1},emdata{1});
set(abshdl,'DisplayName','Absorption','LineWidth',2);
set(emhdl(1),'DisplayName',['Emission, \lambda_{ex} = ',num2str(emwlval(1)),' nm'],'LineWidth',2,'Color',color(1));
axes(axhdl(2));
hold on
%%

for i = 2:length(emwlval)
    emhdl(i)=plot(emwl{i},emdata{i},'LineStyle','-','Color',color(mod(i-1,5)+1),'LineWidth',2,...
        'DisplayName',['Emission, \lambda_{ex} = ',num2str(emwlval(i)),' nm']);
end

%%

for i = 1+length(emwlval):length(emwlval)+length(exwlval)
    emhdl(i)=plot(emwl{i},emdata{i},'LineStyle',':','Color',color(mod(i-1,5)+1),'LineWidth',2,...
        'DisplayName',['Excitation, \lambda_{em} = ',num2str(exwlval(i-length(emwlval))),' nm']);
end

%%
legend('show')
box(axhdl(1),'off')
set(axhdl(1),'YLim',[0 1.2],'YTickMode','auto')
set(axhdl,'XLim',[210 600],'XTickMode','auto')
set(axhdl(2),'YLimMode','auto','YTickMode','auto')

xlabel('Wavelength (nm)')
axes(axhdl(1))
ylabel('Absorbance')
axes(axhdl(2))
ylabel('Emission (cps)')
title(plottitle)
return
%% Get max between two points
[gx,~]=ginput(2);
disp('Select the plot for which you want to get the max');
keyboard;
thisx = get(gco,'XData');
thisy = get(gco,'YData');
[~,~,gli] = closest_member(gx(1),thisx);
[~,~,gri] = closest_member(gx(2),thisx);
gli = gli(1);
gri = gri(1);
if (gli<gri)
    [~,midx]=max(thisy(gli:gri));
    disp(['Max wavelength is ',num2str(thisx(gli+midx-1)),' nm at y = ',num2str(thisy(gli+midx-1))]);
else
    [~,midx]=max(thisy(gri:gli));
    disp(['Max wavelength is ',num2str(thisx(gri+midx-1)),' nm at y = ',num2str(thisy(gri+midx-1))]);
end