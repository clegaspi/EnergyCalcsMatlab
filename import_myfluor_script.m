path = 'C:\Users\clegaspi\Google Drive\Experimental Data\Fluorometer\112113\2TD1X';
fn = '2tdeb';
fignum = 10;
wavenum = 0;

legend_str = 'Ex ';
wl_range = 200:10:450;
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
if (ishandle(fignum))
    clf(fignum,'reset');
end
figure(fignum);
hold on
color = 'rgbkmc';
pattern = {'-',':','--','-.','-',':','--','-.'};


for i = 1:n
    if (wavenum)
        x = 1e7./wl{i};
    else
        x = wl{i};
    end
    
    nidx = find(x == 350);
    x = x(nidx:end);
    y = data{i}(nidx:end);
    
    nval = 1;
    % nidx = find(x == 439);
    % nval = data{i}(nidx);
    nval = max(y);
    
%     if (isempty(nval))
%         continue
%     end
    y=y./nval;
    
    plot(x,y, 'LineStyle', pattern{floor((i-1)/length(color))+1} , 'LineWidth', 2, ...
        'Color',color(mod(i-1,length(color))+1),...
        'DisplayName', [legend_str, num2str(wl_range(i))]);
%    line(repmat(1e7/(1e7/wl_range(i) + 3000),[1 2]),[min(data{i}) max(data{i})]);
    if (i==1)
        legend('show','-DynamicLegend');
    end
end

if (wavenum)
    xlabel('Wavenumbers (cm^{-1})');
else
    xlabel('Wavelength (nm)');
end
ylabel('Intensity (cps)');