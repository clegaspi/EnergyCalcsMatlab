%% Import data
fignum = 11;
wavenum = 0;
if (~exist('abs_last_path', 'var'))
    abs_last_path = pwd;
end
[fn,path,~] = uigetfile([abs_last_path,'\*.csv'],'Select data','MultiSelect','on');

if (~iscell(fn))
    if (fn==0)
        return;
    end
    fn = {fn};
end

abs_last_path = path;

abswl = cell(1,length(fn));
absdata = abswl;

for i = 1:length(fn);
    [abswl{i}, absdata{i}, sampname] = import_uvvis([path,fn{i}]);
end

button = questdlg('Do you need to subtract the blank?','Blank?','No');

if (strcmpi(button,'Yes'))
    [bfn,bpath,~] = uigetfile([abs_last_path,'\*.csv'],'Select data','MultiSelect','off');
    if (bfn==0)
        return;
    end
    
    absbwl = zeros(1,length(bfn));
    absbdata = absbwl;

    [absbwl, absbdata] = import_uvvis([bpath,bfn]);
    
    for i = 1:length(fn);
        if ((max(absbwl)+1e-4 > max(abswl{i})) && (min(absbwl)-1e-4 < min(abswl{i})))
            bldata = interp1(absbwl,absbdata,abswl{i});
            absdata{i} = absdata{i} - bldata;
        else
            throw(MException('UVVis:BadBlankBounds','The wavelength bounds of the blank are not wide enough for this sample!'));
        end
    end    
elseif (strcmpi(button,'Cancel'))
    return;
end





%% Plot 
figure(fignum);
hold on
color = 'rgbkmc';
pattern = {'-',':','--','-.','-',':','--','-.'};

for i = 1:length(fn)
    if (wavenum)
        x = 1e7./abswl{i};
    else
        x = abswl{i};
    end
    plot(x,absdata{i}, 'LineStyle', pattern{floor((i-1)/length(color))+1} , 'LineWidth', 2, ...
        'Color',color(mod(i-1,length(color))+1),...
        'DisplayName', sampname);
end

if (wavenum)
    xlabel('Wavenumbers (cm^{-1})');
else
    xlabel('Wavelength (nm)');
end
ylabel('Absorbance (AU)');
legend('show');