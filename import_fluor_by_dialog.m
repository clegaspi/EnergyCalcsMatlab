%% Import data
fignum = 11;
wavenum = 1;
if (~exist('em_last_path', 'var'))
    em_last_path = pwd;
end
[fn,path,~] = uigetfile([em_last_path,'\*.csv'],'Select data','MultiSelect','on');

if (~iscell(fn))
    if (fn==0)
        return;
    end
    fn = {fn};
end

em_last_path = path;

emwl = cell(1,length(fn));
emdata = emwl;

for i = 1:length(fn);
    [emwl{i}, emdata{i}] = import_fluor([path,fn{i}]);
end

%% Plot 
figure(fignum);
hold on
color = 'rgbkmc';
pattern = {'-',':','--','-.','-',':','--','-.'};

for i = 1:length(fn)
    if (wavenum)
        x = 1e7./emwl{i};
    else
        x = emwl{i};
    end
    plot(x,emdata{i}, 'LineStyle', pattern{floor((i-1)/length(color))+1} , 'LineWidth', 2, ...
        'Color',color(mod(i-1,length(color))+1),...
        'DisplayName', fn{i});
end

if (wavenum)
    xlabel('Wavenumbers (cm^{-1})');
else
    xlabel('Wavelength (nm)');
end
ylabel('Intensity (cps)');
legend('show');