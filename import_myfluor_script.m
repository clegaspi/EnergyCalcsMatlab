uiippath = 'C:\Users\clegaspi\Google Drive\Experimental Data\Fluorometer\092413';
fn = 'methft';
fignum = 12;
wavenum = 0;

legend_str = 'Ex ';
wl_range = 200:10:600;
n = length(wl_range);

%%

if (path(end) ~= '\')
    path = [path, '\'];
end

wl = cell(1,n);
data = cell(1,n);
for i = 1:n
    [wl{i},data{i}]=import_fluor([path,fn,num2str(i-1,'%02u'),'.csv'],2);
    data{i} = smooth(data{i});
%     data{i} = data{i} ./ max(data{i}(44:54));
end

%%
figure(fignum);
hold on
color = 'rgbkmc';
pattern = {'-',':','--','-.','-',':','--','-.'};


for i = 1:10
    if (wavenum)
        x = 1e7./wl{i};
    else
        x = wl{i};
    end
    plot(x,data{i}, 'LineStyle', pattern{floor((i-1)/length(color))+1} , 'LineWidth', 2, ...
        'Color',color(mod(i-1,length(color))+1),...
        'DisplayName', [legend_str, num2str(wl_range(i))]);
    if (i==21)
        legend('show','-DynamicLegend');
    end
end

if (wavenum)
    xlabel('Wavenumbers (cm^{-1})');
else
    xlabel('Wavelength (nm)');
end
ylabel('Intensity (cps)');