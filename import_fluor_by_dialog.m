%% Import data
append_data = 1;
fignum = 5;
wavenum = 0;
if (~exist('em_last_path', 'var'))
    em_last_path = pwd;
end

[tfn,path,~] = uigetfile([em_last_path,'\*.csv'],'Select data','MultiSelect','on');

if (~iscell(tfn))
    if (tfn==0)
        return;
    end
    tfn = {tfn};
end

if (append_data)
    fn = [fn tfn];
else
    fn = tfn;
end

em_last_path = path;

if (~append_data)
    emwl = cell(1,length(fn));
    emdata = emwl;
    start = 1;
else
    start = length(emwl)+1;
    emwl = [emwl cell(1,length(tfn))];
    emdata = [emdata cell(1,length(tfn))];
end

for i = start:(start+length(tfn)-1)
    [emwl{i}, emdata{i}] = import_fluor([path,fn{i}]);
    if (size(emdata{i},2) > 1)
        emdata{i} = mean(emdata{i},2);  % Average if we have multiple runs
        emwl{i} = emwl{i}(:,1); % Collapse wl to 1D
    end
end

%% Plot 
if (~append_data && ishandle(fignum))
    clf(fignum);
    sidx = 1;
    eidx = length(emwl);
elseif (append_data && ishandle(fignum))
    sidx = start;
    eidx = start+length(fn)-1;
else
    sidx = 1;
    eidx = length(emwl);
end

figure(fignum);
hold on
color = 'rgbkmc';
pattern = {'-',':','--','-.','-',':','--','-.'};

for i = sidx:eidx
    if (wavenum)
        x = 1e7./emwl{i};
    else
        x = emwl{i};
    end
    htmp = plot(x,emdata{i}, 'LineStyle', pattern{floor((i-1)/length(color))+1} , 'LineWidth', 2, ...
        'Color',color(mod(i-1,length(color))+1));
    if (length(fn)>=i-sidx+1)
        set(htmp,'DisplayName', fn{i-sidx+1});
    end
end

if (wavenum)
    xlabel('Wavenumbers (cm^{-1})');
else
    xlabel('Wavelength (nm)');
end
ylabel('Intensity (cps)');
if (~append_data)
    legend('show','-DynamicLegend');
end