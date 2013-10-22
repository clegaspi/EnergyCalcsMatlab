% %% DON'T RUN FROM HERE
% 
% while false

%% Run multiple molecules

mols = {'MeLPPP-13mer',.064,[1 186]};
nstates = 25;
norbs = 1000;
poolsize = 7;

run_energy_calcs = false;
load_data = true;
run_struct_opt = false;

if (run_energy_calcs && matlabpool('size') == 0)
    matlabpool(poolsize);
end

for imol = 1:size(mols, 1);
%     system('subst s: /d');
%     system('subst o: /d');
%     system('subst p: /d');
%     
%     rootdir = ['C:\Users\clegaspi\Documents\MATLAB\data\',mols{imol,1},'\'];
%     SFolder = 'Exp';
%     OFolder = 'INDOLib';
%     PFolder = 'GSLib';
% 
%     system(['subst s: "',rootdir,SFolder,'"']);
%     system(['subst o: "',rootdir,OFolder,'"']);
%     system(['subst p: "',rootdir,PFolder,'"']);
    
     field = 0:0.001:mols{imol,2};
     % field = [0:0.001:0.045, 0.0452:0.0002:0.065, 0.066:0.001:mols{imol,2}];
        

%% Are you loading pre-existing data? Set this to false.
if (run_energy_calcs)    
%% 

% There are no PHI angles in this script. To systematically modify the
% dihedral, use BOScriptMultipleAngles.m
    
% waithdl = waitbar(0, [mols{imol,1},' Calculations running...']);
% drawnow;

no_field_dm = [];

if (exist(['j:\NoAngle\',mols{imol,1},'-0VA.mat'], 'file'))
    S = load(['j:\NoAngle\',mols{imol,1},'-0VA.mat']);
    myexp = S.obj;
    if (exist(['k:\', myexp.data(1).indo_hash, '-dm.bin'],'file'))
        no_field_dm = ['k:\', myexp.data(1).indo_hash, '-dm.bin'];
    end
elseif (dblcmp(field(1),0))
    k = 1;
    disp(num2str(field(k)));

    myece = ECEParams('AM1', {}, true, field(k), norbs, nstates,...
        [pwd,'\data\params_for_all.txt']);

    myexp = EnergyCalcExp(myece,...
        ['j:\NoAngle\',mols{imol,1},'.dat'],...
        mols{imol,3},...
        ['j:\NoAngle\',mols{imol,1},'-',num2str(field(k)),'VA.mat'],...
        'k:\',...
        'l:\',...
        false,...
        ['MeLPPP-13mer Field calcs - ',num2str(field(k)),'VA']);
    myexp.run('quiet');
    
    no_field_dm = ['k:\', myexp.data(1).indo_hash, '-dm.bin'];
end

parfor pk = 2:numel(field)
    disp(num2str(field(pk)));
%     waithdldata = get(waithdl,'userdata');
%     if (waithdldata.bail == 0)
%         waitbar(n / (numel(phi1)*numel(phi2)*numel(field)), waithdl);
%     else
%         keyboard;
%     end

    if (exist(['j:\NoAngle\',mols{imol,1},'-',num2str(field(pk)),'VA.mat'], 'file'))
        continue;
    else
        myece = ECEParams('AM1', {}, true, field(pk), norbs, nstates,...
            [pwd,'\data\params_for_all.txt']);

        if (~isempty(no_field_dm))
            myece.dm_guess = no_field_dm;
            myece.try_default_first = true;
        end

        pmyexp = EnergyCalcExp(myece,...
            ['j:\NoAngle\',mols{imol,1},'.dat'],...
            mols{imol,3},...
            ['j:\NoAngle\',mols{imol,1},'-',num2str(field(pk)),'VA.mat'],...
            'k:\',...
            'l:\',...
            false,...
            ['MeLPPP-13mer Field Calcs - ',num2str(field(pk)),'VA']);
        pmyexp.run('quiet');
    end
%     if (~exist('n','var'))
%         n=0;
%     end
%     n = n + 1;
end

%% Doing calcs? Terminate
end
% close(waithdl,'force');
%% Skip data loading because you're doing BOBL opt? Set to false
if (load_data)
    
%% Load data that's been calculated already

dp = zeros(1,length(field));
dpgs = zeros(1,length(field));
dp2exc = zeros(1,length(field));
dp3exc = zeros(1,length(field));
dp4exc = zeros(1,length(field));
dp5exc = zeros(1,length(field));
dpmag = zeros(1,length(field));
hlgap = zeros(1,length(field));
% eexc = zeros(length(field),nstates);
eexc = zeros(1,length(field)*nstates);
igs = zeros(1,length(field));
opint = zeros(length(field),nstates);
tpint = zeros(length(field),nstates);

exdp = zeros(1,length(field));
exdpgs = zeros(1,length(field));
exdp2exc = zeros(1,length(field));
exhlgap = zeros(1,length(field));
exeexc = zeros(length(field),nstates);
exigs = zeros(1,length(field));
exdeltaAM1 = zeros(1,length(field));
exopint = zeros(length(field),nstates);
extpint = zeros(length(field),nstates);

waithdl = waitbar(0, [mols{imol,1},' Loading data...']);
drawnow;

for idx = 1:numel(field)
    waithdldata = get(waithdl,'userdata');
    if (waithdldata.bail == 0)
        waitbar((idx-1)/numel(field), waithdl);
        drawnow;
    else
        keyboard;
    end

    S = load(['j:\NoAngle\',mols{imol,1},'-',num2str(field(idx)),'VA.mat']);
    myexp = S.obj;
    
    myexp.data(1).update_paths('l:\','k:\');
    
    tmp = myexp.get_field('indo.dipole',1,1);
    dpgs(idx) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = myexp.get_field('indo.dipole',2,2);
    dp(idx) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = myexp.get_field('indo.dipole',3,3);
    dp2exc(idx) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = myexp.get_field('indo.dipole',4,4);
    dp3exc(idx) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = myexp.get_field('indo.dipole',5,5);
    dp4exc(idx) = sum(tmp .^ 2, 1) .^ (0.5);
    
    tmp = myexp.get_field('indo.dipole',6,6);
    dp5exc(idx) = sum(tmp .^ 2, 1) .^ (0.5);

%     maxtpstate = 16;    % This will mean nothing until we see the data
%     tmp = myexp.get_field('indo.dipole',maxtpstate,maxtpstate);
%     dpmag(idx) = sum(tmp .^ 2, 1) .^ (0.5);

    nfill = myexp.get_field('indo.nfilled');
    tmp = myexp.get_field('indo.orbE',[nfill nfill+1]);
    hlgap(idx) = tmp(2) - tmp(1);
    
    eexc(idx,1:nstates) = myexp.get_field('Eexc',:);
    igs(idx) = myexp.get_field('indo.esci',1);
    opint(idx,1:nstates) = myexp.get_field('Tint',:);
%     tmp = squeeze(myexp.get_field('indo.r',2,:,:));  % From n=2
%     midx = 2;  % for n=2
%     for l = 1:nstates
%         tmpeexc = eexc(idx,1:nstates);
%         [~,midx] = max(opint(idx,1:nstates));
%         tmp = squeeze(myexp.get_field('indo.r',midx,:,:));  % From intense OP
%         tpint(idx,l) = (tmpeexc(l) - tmpeexc(midx)) * sum(tmp(l,:) .^ 2);
%     end
%     
end
close(waithdl,'force');

% Make all cells because cross-compatibility with code
% Ground state structure
dp = {dp};
dpgs = {dpgs};
hlgap = {hlgap};
eexc = {eexc};
igs = {igs};
opint = {opint};
tpint = {tpint};
dp2exc = {dp2exc};
dp3exc = {dp3exc};
dp4exc = {dp4exc};
dp5exc = {dp5exc};
dpmag = {dpmag};

% Excited state structure
exdp = {exdp};
exdpgs = {exdpgs};
exhlgap = {exhlgap};
exeexc = {exeexc};
exigs = {exigs};
exdeltaAM1 = {exdeltaAM1};
exopint = {exopint};
extpint = {extpint};
exdp2exc = {exdp2exc};


end

end

%% Skip structure opt? Debugging?

if (run_struct_opt)

%% Do structure optimizations

% fields_to_opt_struct = field;
fields_to_opt_struct = 0;

boandrdata = cell(1,1);

for fstrn = fields_to_opt_struct;
    k = find(arrayfun(@(x)dblcmp(x,fstrn), field));
    
    S = load(['j:\NoAngle\',mols{imol,1},'-',num2str(fstrn),'VA.mat']);
    myexp = S.obj;
    
    thisdat = myexp.data(1);

    thisdat.load_to_memory('indo','load');
    indo = thisdat.raw_indo;
    efield = indo.config.field;

    thisdat.generate_ampac_file('out', 'l:\');

    [rgs, rex, deltar, eszmat, ampac_energy, exindo, ~, numruns] = OptExcStateStructure('l:\', thisdat.ampac_hash,...
        'k:\', thisdat.indo_hash,...
        'field', efield,...
        'indo', indo,...
        'algorithm','paulingwitham1',...
        'nstates', nstates);

    % Second excited state
%                 [rgs, rex, deltar, ~, numruns, eszmat, ampac_energy] = OptExcStateStructure('l:\', myexp.data(1).ampac_hash,...
%                     'k:\', myexp.data(1).indo_hash,...
%                     'field', efield,...
%                     'indo', indo,...
%                     'state', 3,...
%                     'readifexist');

%                 if (isempty(boandrdata))
%                     natom = myexp.get_field('ampac.natom');
%                     boandrdata = cell(natom,6);
%                 end

    boandrdata{1,1}{k, 1} = fstrn;
    boandrdata{1,1}{k, 2} = rgs;
    boandrdata{1,1}{k, 3} = rex;
    boandrdata{1,1}{k, 4} = deltar;
    boandrdata{1,1}{k, 5} = numruns;
    boandrdata{1,1}{k, 6} = eszmat;
    % idx = idx + 1;

    exdeltaAM1{1,1}(k) = ampac_energy - myexp.data(1).Ehf;

    myexp.data(1).load_to_memory('ampac','load');
    amp = myexp.data(1).raw_ampac;
    amp.ampac_succeed = true;

    exindo = Indo.LoadExistingData(['k:\', myexp.data(1).indo_hash, '-new.ido'],[],[],[]);

    tmp = exindo.dipole(1,1);
    exdpgs{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = exindo.dipole(2,2);
    exdp{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = exindo.dipole(3,3);
    exdp2exc{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    exnfill = exindo.nfilled;
    tmp = exindo.orbE([exnfill exnfill+1]);
    exhlgap{1,1}(k) = tmp(2) - tmp(1);

    tempds = ECEDataStruct(amp, exindo,[], [], [], [], []);
    exeexc{1,1}(k,1:nstates) = tempds.Eexc(:);

    exigs{1,1}(k) = exindo.esci(1);

    exopint{1,1}(k,1:nstates) = tempds.Tint(:);
    tmp = squeeze(exindo.r(2,:,:));
    for l = 1:nstates
        extpint{1,1}(k,l) = (exeexc{1,1}(k,l) - exeexc{1,1}(k,2)) * sum(tmp(l,:) .^ 2);
    end
end
    


%% Save calcs for multiple molecules

save(['j:\',mols{imol,1},'-BOScript-PaulingwithAM1.mat']);

end

disp('Calculation finished');
return;

%% GS Structure Energy vs. Field
figure(8)
hold on

xaxis = field;
% xaxis = 0;

if (~exist('nstates','var'))
    nstates = length(eexc{1,1}(k,:));
end

for k = 1:length(xaxis)
    plot(xaxis(k), eexc{1,1}(k,1:nstates) + igs{1,1}(k) - igs{1,1}(1), 'g^');
end

maxopint = max(opint{1,1}(:));
for k = 1:length(xaxis)
    for l = 1:nstates
        if (opint{1,1}(k,l)*30/maxopint > 1)
            plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k) - igs{1,1}(1), 'bo', 'MarkerSize', (opint{1,1}(k,l)*30 / maxopint) + 1e-3);
        end
    end
end

xlabel('Applied Field Strength / [V cm^-^1] x 10^-^8');
ylabel('Relative Energy / [eV]')

% maxtpint = max(abs(tpint{1,1}(:)));
% 
% for k = 1:length(xaxis)
%     for l = 1:nstates
%         if (tpint{1,1}(k,l) > 0)
%             plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k) - igs{1,1}(1), 'rs', 'MarkerSize', (tpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
%         elseif (tpint{1,1}(k,l) < 0)
%             plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k) - igs{1,1}(1), 'rd', 'MarkerSize', (-tpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
%         end
%     end
% end

%% Draw grid lines on points

figure(8);
hold on;

grid_hdl = zeros(1,length(xaxis));
for i = 1:length(xaxis)
    grid_hdl(i) = plot([xaxis(i) xaxis(i)], [(eexc{1,1}(i,1) + igs{1,1}(i)) (eexc{1,1}(i,nstates) + igs{1,1}(i))] - igs{1,1}(1),...
        'k:');
end

%% Draw line to follow states only, calculate dipole moment from slope

if (exist('lnhnd','var') && ishandle(lnhnd))
    delete(lnhnd);
end

figure(8);
hold on;


% [x,y] = ginput(2);
% m = (y(2) - y(1)) / (x(2) - x(1));
% b = y(1) - m*x(1);

[x,y] = ginput(3);
fit = polyfit(x,y,1);
m = fit(1);
b = fit(2);

lnhnd = plot(xaxis, m .* xaxis + b, 'k-');

disp(['The equation of the line is ',num2str(m),'x + ',num2str(b)]);
disp(['The dipole moment of this state is ',num2str(-m),' e-A (',num2str(-m*4.802456),' D)']);

%% Draw two lines to follow two states and calculate their intersection

if (exist('lnhnd1','var') && ishandle(lnhnd1))
    delete(lnhnd1);
end

if (exist('lnhnd2','var') && ishandle(lnhnd2))
    delete(lnhnd2);
end


figure(8);
hold on;

% [x,y] = ginput(4);
% m1 = (y(2) - y(1)) / (x(2) - x(1));
% b1 = y(1) - m1*x(1);
% 
% m2 = (y(4) - y(3)) / (x(4) - x(3));
% b2 = y(3) - m2*x(3);

[x,y] = ginput(6);
fit = polyfit(x(1:3),y(1:3),1);
m1 = fit(1);
b1 = fit(2);
fit = polyfit(x(4:6),y(4:6),1);
m2 = fit(1);
b2 = fit(2);

lnhnd1 = plot(xaxis, m1 .* xaxis + b1, 'c-');
lnhnd2 = plot(xaxis, m2 .* xaxis + b2, 'm-');

interx = (b2-b1)/(m1-m2);
intery = m1*interx + b1;

disp(['The intersection is ( ',num2str(interx),' , ',num2str(intery),' )']);


%% Remove lines
delete(lnhnd1,lnhnd2);

%% Get slope and y-int of nBu. Calculate appx value of exciton binding energy
figure(8);  % Need to have energy vs. field plotted here already

get_slope_from_n2_after_cross = true;
use_n2_at_zero_for_eb = true;
hold_nofield_energy = false;

if (exist('ebhnd','var') && ishandle(ebhnd))
    delete(ebhnd);
end

figure(8);
hold on;

if (get_slope_from_n2_after_cross)
    dcm_obj = datacursormode(8);
    for i = 1:2
        disp(['Pick point ',num2str(i),' with the data cursor and then type return']);
        figure(8);
        keyboard;
        info_struct = getCursorInfo(dcm_obj);
        pts(i) = info_struct.Position(1);
    end
    
    start_idx = find(field == pts(1));
    end_idx = find(field == pts(2));
    
    x = field(start_idx:end_idx)';
    y = eexc{1}(start_idx:end_idx,2) + igs{1,1}(start_idx:end_idx)' - igs{1,1}(1);
else
    [x,y] = ginput(3);
end

fit = polyfit(x,y,1);
m = fit(1);
b = fit(2);

if (~hold_nofield_energy || ~exist('eby','var'))
    if (use_n2_at_zero_for_eb)
        eby = eexc{1}(1,2);
    else
        [~,eby] = ginput(1);
    end
end

ebenergy = b-eby;

ebhnd = plot(xaxis, m .* xaxis + b, 'k-');

disp(['The equation of the line is ',num2str(m),'x + ',num2str(b)]);
disp(['The dipole moment of this state is ',num2str(-m),' e-A (',num2str(-m*4.802456),' D)']);
disp(['The exciton binding energy is approximately ', num2str(ebenergy), ' eV']);

%% Plot ES Energy Data

figure(8)
hold on

% xaxis = 0:0.01:0.05;
i = 1; j = 1;
%xaxis = [0:0.1:(floor(10*maxdata(p1i(i),p2i(j)) - 0.1) / 10) secondthreshold(p1i(i),p2i(j)):0.1:(secondthreshold(p1i(i),p2i(j))+1)];
% xaxis = xaxis(1:8);
xaxis = field;

for k = 1:length(xaxis)
    plot(xaxis(k), exeexc{1,1}(k,:) + igs{1,1}(k) + exdeltaAM1{1,1}(k), 'k^');
end

% maxexopint = max(exopint{1,1}(:));
for k = 1:length(xaxis)
    for l = 1:nstates
        plot(xaxis(k), exeexc{1,1}(k,l) + igs{1,1}(k) + exdeltaAM1{1,1}(k), 'mo', 'MarkerSize', (exopint{1,1}(k,l)*30 / maxopint) + 1e-3);
    end
end

% maxtpint = max(tpint{1,1}(:));

for k = 1:length(xaxis)
    for l = 1:nstates
        if (tpint{1,1}(k,l) > 0)
            plot(xaxis(k), exeexc{1,1}(k,l) + igs{1,1}(k) + exdeltaAM1{1,1}(k), 'cs', 'MarkerSize', (extpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
        end
    end
end

% %% BEGIN RUNNING HERE
% end
%% Draw structure with delta r color coding

field_strength = 0.101;
mol_name = '8-merPPV';

S = load(['j:\NoAngle\',mol_name,'-0VA.mat']);
myexp = S.obj;

myexp.data(1).generate_ampac_file('out');

[xyz, atom_type] = ampac_to_xyz(['l:\', myexp.data(1).ampac_hash, '.out']);
[exxyz, ~] = ampac_to_xyz(boandrdata{1,1}{find(arrayfun(@(x)dblcmp(x, field_strength), field)), 6}); 
atomsbonds = get_connectivity(['l:\', myexp.data(1).ampac_hash, '.out']);

cidx = find(cellfun(@(x)strcmpi(x,'c'), atom_type));

% xyz = xyz(cidx,:);
% exxyz = exxyz(cidx,:);

figure(7);
hold on;
axis equal;
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) min(xyz(:,3))-1 max(xyz(:,3))+1]);
camorbit(0,270);
for i = 1:size(xyz,1)
    scatter3(xyz(i,1),xyz(i,2),xyz(i,3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k');
end

dr = zeros(1,length(atomsbonds));
dridx = 1;

for i = cidx'
    btd = intersect(atomsbonds{i}, cidx);
    for j = btd 
        dat = [xyz(i,:)' xyz(j,:)'];
        exdat = [exxyz(i,:)' exxyz(j,:)'];
        dr(dridx) = (sum((exdat(:,1)-exdat(:,2)) .^ 2) ^ 0.5) - (sum((dat(:,1)-dat(:,2)) .^ 2) ^ 0.5);
        dridx = dridx + 1;
    end
end

dridx = 1;
maxdr = max(abs(dr));

for i = cidx'
    btd = intersect(atomsbonds{i}, cidx);
    for j = btd 
        dat = [xyz(i,:)' xyz(j,:)'];
        
        if (dr(dridx) < 0)
            lnclr = 'r';
        elseif (dr(dridx) > 0)
            lnclr = 'c';
        else
            lnclr = 'b';
        end
        
        mag = (abs(dr(dridx)) / maxdr) * 7 + 1e-4;
        
        line(dat(1,:),dat(2,:),dat(3,:),'color',lnclr,'linewidth',mag)
        dridx = dridx + 1;
    end
end

%% Dipole moments
figure(17);
hold on;
xaxis = field;
k = 1;
scatter(xaxis,reshape(dpgs{k}(1:length(xaxis)) .* 4.8,1,[]),'MarkerEdgeColor','g','Marker','*');
scatter(xaxis,reshape(dp{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','b','Marker','o');
scatter(xaxis,reshape(dp2exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','r','Marker','+');
% scatter(xaxis,reshape(dp3exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','k','Marker','d');
% scatter(xaxis,reshape(dp4exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','c','Marker','^');
% scatter(xaxis,reshape(dp5exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','m','Marker','x');

xlabel('Applied Field / [V cm^-^1] x 10^-^8');
ylabel('Dipole Moment / [D]')
legend('n = 1','n = 2','n = 3','n = 4','n = 5','n = 6');

%% Bond Length Alternation Plot
% Enter atom numbers of S-D-S units where 1-2 and 3-4 are singles and 2-3
% is a double
bla_anums = [3 7 8 9; 14:17; 22:25; 30:33; 38:41; 46:49; 54:57];    % 8-merPPV
field_strength = 0.101;

[xyz, atom_type] = ampac_to_xyz(['l:\', myexp.data(1).ampac_hash, '.out']);
[exxyz, ~] = ampac_to_xyz(boandrdata{1,1}{arrayfun(@(x)dblcmp(x, field_strength), field), 6}); 

gs_bla = zeros(size(bla_anums, 1),1);
ex_bla = gs_bla;

for i = 1:length(gs_bla)
    sb1 = sum((xyz(bla_anums(i,1),:) - xyz(bla_anums(i,2),:)) .^ 2) ^ 0.5;
    sb2 = sum((xyz(bla_anums(i,3),:) - xyz(bla_anums(i,4),:)) .^ 2) ^ 0.5;
    db = sum((xyz(bla_anums(i,2),:) - xyz(bla_anums(i,3),:)) .^ 2) ^ 0.5;
    gs_bla(i) = ((sb1 + sb2) / 2) - db;
    
    sb1 = sum((exxyz(bla_anums(i,1),:) - exxyz(bla_anums(i,2),:)) .^ 2) ^ 0.5;
    sb2 = sum((exxyz(bla_anums(i,3),:) - exxyz(bla_anums(i,4),:)) .^ 2) ^ 0.5;
    db = sum((exxyz(bla_anums(i,2),:) - exxyz(bla_anums(i,3),:)) .^ 2) ^ 0.5;
    ex_bla(i) = ((sb1 + sb2) / 2) - db;
end

figure(5);
hold on;
% plot(1:length(gs_bla), gs_bla, '-g');
plot(1:length(ex_bla), ex_bla, ':r');
xlabel('Repeat Unit (n)');
ylabel('Bond Length Alternation (Angstroms)');

%% Charge density plot

fignum = 5;
exciton = false;
mol_name = 'MeLPPP-13mer';

pick_a_point = false;
% If pick_a_point, the following will be overwritten
field_strength = 0;
state_n = 2;

% Put how your y-axis is plotted here. Otherwise set y_offset = 0 for
% absolute energies.
S = load(['j:\NoAngle\',mol_name,'-0VA.mat']);
myexp = S.obj;
myexp.data(1).update_paths('l:\','k:\');
myexp.data(1).load_to_memory('indo','load');
myindo = myexp.data(1).raw_indo;
y_offset = myindo.esci(1);

centroid_offset = [0 0 0];

% Graphing charge density
binsize = 2;    % Angstroms. +- amount for bin capture
dr = 0.05;  % How many angstroms to incrememnt for bin center
nsmooth = 10;  % How many times to smooth the curve
triangle = true;    % Triangle (weighted) smooth or rectangular (average)
tsize = 15;  % How big is the triangle?
baraswell = false;   % Plot the bar graph on top of the line plot
plot_line_color = 'b';

if (pick_a_point)
    dcm_obj = datacursormode(8);
    info_struct = getCursorInfo(dcm_obj);
    while (isempty(info_struct))
        disp('Pick a point with the data cursor and then type return');
        figure(8);
        keyboard;
        info_struct = getCursorInfo(dcm_obj);
    end        
    field_strength = info_struct.Position(1);
end


S = load(['j:\NoAngle\',mol_name,'-',num2str(field_strength),'VA.mat']);
myexp = S.obj;
myexp.data(1).update_paths('l:\','k:\');
myexp.data(1).load_to_memory('indo','load');
myindo = myexp.data(1).raw_indo;

if (pick_a_point)
    found_it = false;
    for state_n = 1:length(myindo.esci)
        if (dblcmp(myindo.esci(state_n) - y_offset, info_struct.Position(2), 1e-4))
            found_it = true;
            break;
        end
    end
end

if (pick_a_point && ~found_it)
    throw(MException('ChgDensity:NoPointFound','Charge Density Plotting: No point found to plot!'));
end

disp(['Plotting ( ', num2str(field_strength),' , ',num2str(myindo.esci(state_n)-y_offset),' ) (n=',num2str(state_n),')']);

% indo_filepath = ['k:\',myexp.data(1).indo_hash,'-dm.bin'];
% fid = fopen(indo_filepath);
% matdim = fread(fid,1,'int');
% gsdm = fread(fid,[matdim,matdim],'double');
% fclose(fid);

if (exciton)
    deltadm = diffDensity(myindo, state_n, 'exciton');    % Look at exciton
else
    deltadm = diffDensity(myindo, state_n);     % Look at charge probability
end

% deltadm = gsdm;

% esdm = gsdm + deltadm;

% Get cartesian coordinates from ampac
myexp.data(1).generate_ampac_file('out','l:\');
[xyz, atom_type] = ampac_to_xyz(['l:\', myexp.data(1).ampac_hash, '.out']); % Ground state
zmat = ampac_to_zmatrix(['l:\', myexp.data(1).ampac_hash, '.out']);

% Mulliken charge
mchg = zeros(1,length(xyz));

for i = 1:length(myindo.aorbAtom)
    mchg(myindo.aorbAtom(i)) = mchg(myindo.aorbAtom(i)) - deltadm(i,i);
end

for i = 1:size(zmat,1)
    if (strcmpi(zmat{i,2},'h'))
        carbon_att = str2double(zmat{i,6});
        % mchg(carbon_att) = mchg(carbon_att) + 1 - mchg(i);
        mchg(carbon_att) = mchg(carbon_att) + mchg(i);
        mchg(i) = 0;
    end
end

% Rotate the Cartesian coordinates so x axis is long axis
f_atom = myexp.axis_params(1);
e_atom = myexp.axis_params(2);
oldv = xyz(e_atom,:) - xyz(f_atom,:);   % Get vector between front and end atom if front atom is at origin
rotmat = vrrotvec2mat(vrrotvec(oldv,[1 0 0]));  % Rotation matrix to have the long axis along x

rot_xyz = xyz;

centroid = zeros(1,3);  % Calculate the centroid for transition dipole vector arrow

for i = 1:length(xyz(:,1))
    rot_xyz(i,:) = xyz(i,:) - xyz(f_atom,:);    % Translation to move f_atom to origin
    rot_xyz(i,:) = rotmat * rot_xyz(i,:)';      % Rotation to move long axis to x axis
    centroid = centroid + rot_xyz(i,:);     % Centroid is the average of all coordinates
end

centroid = centroid ./ length(xyz(:,1));

if (ishandle(fignum))
    clf(fignum,'reset');
end
if (exist('harr','var') && ishandle(harr))
    delete(harr);
end
figure(fignum)
hold on;
axis equal;
axis([min(rot_xyz(:,1))-1 max(rot_xyz(:,1))+1 min(rot_xyz(:,2))-1 max(rot_xyz(:,2))+1 min(rot_xyz(:,3))-1 max(rot_xyz(:,3))+1]);
cptcmap('scaled_mulliken','mapping','direct')

scatter3(rot_xyz(:,1),rot_xyz(:,2),rot_xyz(:,3),120,mchg*1e3/max(abs(mchg)),'o','filled','MarkerEdgeColor','k')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
xlabel('x');
ylabel('y');
zlabel('z');

% Graph the arrow showing GS transition dipole (arrow size does not correlate
% to magnitude
tdv = reshape(myindo.r(1,state_n,:),1,3);
tdv = tdv .* (norm(centroid)/2) ./ norm(tdv);
if (tdv(2) > 0)
    carr = [0 1 0];
else
    carr = [1 0 0];
end
harr = arrow(centroid+centroid_offset,centroid'+centroid_offset'+rotmat*tdv','Width',3,'FaceColor',carr);

% Graph of spatial distribution of charges

shift_xyz = rot_xyz(:,1)-min(rot_xyz(:,1));
graph_pts = 0:dr:(max(shift_xyz)-binsize);
chg_in_bin = zeros(size(graph_pts));

for i = 1:length(graph_pts)
    for j = 1:length(shift_xyz)
        if (shift_xyz(j) >= (graph_pts(i) - binsize) && ...
                shift_xyz(j) <= (graph_pts(i) + binsize))
            chg_in_bin(i) = chg_in_bin(i) + mchg(j);
        end
    end
end

% Curve smoothing
for idx = 1:nsmooth
    if (triangle)
        for i = (1+tsize):(length(chg_in_bin)-tsize)
            tmp = ((tsize+1)/((tsize+1)^2))*chg_in_bin(i);
            for j = 1:tsize
                tmp = tmp + (j/((tsize+1)^2))*chg_in_bin(i-tsize+j-1) + ...
                    (j/((tsize+1)^2))*chg_in_bin(i+tsize-j+1);
            end
            chg_in_bin(i) = tmp;
        end
        
        % Change tsize for the end points
        thist = tsize;
        for i = tsize:-1:2
            thist = thist - 1;
            tmp = ((thist+1)/((thist+1)^2))*chg_in_bin(i);
            for j = 1:thist
                tmp = tmp + (j/((thist+1)^2))*chg_in_bin(i-thist+j-1) + ...
                    (j/((thist+1)^2))*chg_in_bin(i+thist-j+1);
            end
            chg_in_bin(i) = tmp;
            
            tmp = ((thist+1)/((thist+1)^2))*chg_in_bin(end-i+1);
            for j = 1:thist
                tmp = tmp + (j/((thist+1)^2))*chg_in_bin(end-i+1-thist+j-1) + ...
                    (j/((thist+1)^2))*chg_in_bin(end-i+1+thist-j+1);
            end
            chg_in_bin(end-i+1) = tmp;
        end
    else
        for i = 1:length(chg_in_bin)-1
            chg_in_bin(i) = (chg_in_bin(i) + chg_in_bin(i+1)) / 2;
        end
    end
end

figure(234)
if (exciton)
    chg_in_bin = -chg_in_bin;
end
plot(graph_pts,chg_in_bin,'Color',plot_line_color)
ylabel('Charge Density')
if (baraswell)
    bin_xyz=ceil((rot_xyz(:,1)-min(rot_xyz(:,1)))./binsize + 1);
    
    bar_bin = zeros(max(bin_xyz),1);
    for i=1:length(bin_xyz)
        bar_bin(bin_xyz(i)) = bar_bin(bin_xyz(i)) + mchg(i);
    end
    
    figure(234)
    hold on
    if (exciton)
        bar_bin = -bar_bin;
    end
    bar([0:length(bar_bin)-1]*binsize, bar_bin)
    hold off
    ylabel('Electron Density')
end
xlabel('x coordinate / [m] x 10^-^1^0')

disp(['The dipole moment of this structure is ', num2str(sum(myindo.r(state_n,state_n,:).^2)^0.5),...
    ' e-A (', num2str(4.8*sum(myindo.r(state_n,state_n,:).^2)^0.5), ' D)']);
disp(['The normalized transition dipole vector is (',num2str(myindo.r(1,state_n,:)./norm(reshape(myindo.r(1,state_n,:),1,3))),')']);

    %% Get center-to-center distance of bipolarons

figure(234)
[x,~] = ginput(4);  % Click l&r side of l-h max and of r-h max

for i = 1:4
    [~,~,tmp] = closest_member(x(i), graph_pts);
    x(i) = tmp(1);
end

[~,lhs] = max(chg_in_bin(x(1):x(2)));
[~,rhs] = max(chg_in_bin(x(3):x(4)));

sep = graph_pts(rhs+x(3)-1) - graph_pts(lhs+x(1)-1);

disp(['Separation distance is ', num2str(sep), ' Angstroms']);


