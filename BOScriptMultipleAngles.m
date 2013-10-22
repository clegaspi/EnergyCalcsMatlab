% %% DON'T RUN FROM HERE
% 
% while false

%% Run multiple molecules

mols = {'8-merPPV',0.117,[6 62]};

for imol = 1:size(mols, 1);
%     system('subst s: /d');
%     system('subst o: /d');
%     system('subst p: /d');
%     
%     rootdir = ['F:\ProtonCalcs\',mols{imol,1},'\'];
%     SFolder = 'Exp';
%     OFolder = 'INDOLib';
%     PFolder = 'GSLib';
% 
%     system(['subst s: "',rootdir,SFolder,'"']);
%     system(['subst o: "',rootdir,OFolder,'"']);
%     system(['subst p: "',rootdir,PFolder,'"']);
%     
     field = 0:0.001:mols{imol,2};
        
    
%% 

phi1 = 90:5:180;
phi2 = 180;

boandrdata = cell(length(phi1),length(phi2));

% Ground state structure
dp = cell(length(phi1),length(phi2));
dpgs = cell(length(phi1),length(phi2));
hlgap = cell(length(phi1),length(phi2));
eexc = cell(length(phi1),length(phi2));
igs = cell(length(phi1),length(phi2));
opint = cell(length(phi1),length(phi2));
tpint = cell(length(phi1),length(phi2));
dp2exc = cell(length(phi1),length(phi2));
dp3exc = cell(length(phi1),length(phi2));
dp4exc = cell(length(phi1),length(phi2));
dpmag = cell(length(phi1),length(phi2));

% Excited state structure
exdp = cell(length(phi1),length(phi2));
exdpgs = cell(length(phi1),length(phi2));
exhlgap = cell(length(phi1),length(phi2));
exeexc = cell(length(phi1),length(phi2));
exigs = cell(length(phi1),length(phi2));
exdeltaAM1 = cell(length(phi1),length(phi2));
exopint = cell(length(phi1),length(phi2));
extpint = cell(length(phi1),length(phi2));
exdp2exc = cell(length(phi1),length(phi2));


% waithdl = waitbar(0, [mols{imol,1},' Calculations running...']);

for k = 1:numel(field)
%     waithdldata = get(waithdl,'userdata');
%     if (waithdldata.bail == 0)
%         waitbar((k-1) / (numel(phi1)*numel(phi2)*numel(field)), waithdl);
%     else
%         keyboard;
%     end

    if (exist(['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat'], 'file'))
        S = load(['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat']);
        myexp = S.obj;
    else
        myece = ECEParams('AM1', {phi1}, true, field(k), 800, 25,...
            's:\1Angle\params.txt');

        if (field(k) ~= 0)
            myece.dm_guess = ['o:\', myexp.data(1).indo_hash, '-dm.bin'];
            myece.try_default_first = true;
        end

        myexp = EnergyCalcExp(myece,...
            ['s:\1Angle\',mols{imol,1},'.dat'],...
            mols{imol,3},...
            ['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat'],...
            'o:\',...
            'p:\',...
            false,...
            ['Bond Order calcs - ',num2str(field(k)),'VA']);
        myexp.run('quiet');
    end
    
    for i = 1:numel(phi1)
        for j = 1:numel(phi2)
            if (k == 1)
                dp{i,j} = zeros(1,length(field));
                dpgs{i,j} = zeros(1,length(field));
                dp2exc{i,j} = zeros(1,length(field));
                dp3exc{i,j} = zeros(1,length(field));
                dp4exc{i,j} = zeros(1,length(field));
                dpmag{i,j} = zeros(1,length(field));
                hlgap{i,j} = zeros(1,length(field));
                eexc{i,j} = zeros(length(field),25);
                igs{i,j} = zeros(1,length(field));
                opint{i,j} = zeros(length(field),25);
                tpint{i,j} = zeros(length(field),25);

                exdp{i,j} = zeros(1,length(field));
                exdpgs{i,j} = zeros(1,length(field));
                exdp2exc{i,j} = zeros(1,length(field));
                exhlgap{i,j} = zeros(1,length(field));
                exeexc{i,j} = zeros(length(field),25);
                exigs{i,j} = zeros(1,length(field));
                exdeltaAM1{i,j} = zeros(1,length(field));
                exopint{i,j} = zeros(length(field),25);
                extpint{i,j} = zeros(length(field),25);
            end

            % For PHI1 only
            tmp = myexp.get_field('indo.dipole',1,1,phi1(i));
            dpgs{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',2,2,phi1(i));
            dp{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',3,3,phi1(i));
            dp2exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',4,4,phi1(i));
            dp3exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',5,5,phi1(i));
            dp4exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            maxtpstate = 16;    % This will mean nothing until we see the data
            tmp = myexp.get_field('indo.dipole',maxtpstate,maxtpstate,phi1(i));
            dpmag{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

            nfill = myexp.get_field('indo.nfilled',phi1(i));
            tmp = myexp.get_field('indo.orbE',[nfill nfill+1],phi1(i));
            hlgap{i,j}(k) = tmp(2) - tmp(1);

            eexc{i,j}(k,1:25) = myexp.get_field('Eexc',:,phi1(i));
            igs{i,j}(k) = myexp.get_field('indo.esci',1,phi1(i));
            opint{i,j}(k,1:25) = myexp.get_field('Tint',:,phi1(i));
            tmp = squeeze(myexp.get_field('indo.r',2,:,:,phi1(i)));
            for l = 1:25
                tpint{i,j}(k,l) = (eexc{i,j}(k,l) - eexc{i,j}(k,2)) * sum(tmp(l,:) .^ 2);
            end

            %For PHI1 and PHI2

%             tmp = myexp.get_field('indo.dipole',1,1,phi1(i),phi2(j));
%             dpgs{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%             
%             tmp = myexp.get_field('indo.dipole',2,2,phi1(i),phi2(j));
%             dp{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%             
%             tmp = myexp.get_field('indo.dipole',3,3,phi1(i),phi2(j));
%             dp2exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%             
%             tmp = myexp.get_field('indo.dipole',4,4,phi1(i),phi2(j));
%             dp3exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%         
%             tmp = myexp.get_field('indo.dipole',5,5,phi1(i),phi2(j));
%             dp4exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%             
%             maxtpstate = 16;    % This will mean nothing until we see the data
%             tmp = myexp.get_field('indo.dipole',maxtpstate,maxtpstate,phi1(i),phi2(j));
%             dpmag{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);
%             
%             nfill = myexp.get_field('indo.nfilled',phi1(i),phi2(j));
%             tmp = myexp.get_field('indo.orbE',[nfill nfill+1],phi1(i),phi2(j));
%             hlgap{i,j}(k) = tmp(2) - tmp(1);
%             
%             eexc{i,j}(k,1:25) = myexp.get_field('Eexc',:,phi1(i),phi2(j));
%             igs{i,j}(k) = myexp.get_field('indo.esci',1,phi1(i),phi2(j));
%             opint{i,j}(k,1:25) = myexp.get_field('Tint',:,phi1(i),phi2(j));
%             tmp = squeeze(myexp.get_field('indo.r',2,:,:,phi1(i),phi2(j)));
%             for l = 1:25
%                 tpint{i,j}(k,l) = (eexctemp(k,l) - eexctemp(k,2)) * sum(tmp(l,:) .^ 2);
%             end
        end
    end
end

% close(waithdl,'force');

%% Run structure calcs
% PHI angle values and fields to opt structure {phi1,phi2,fields}
phifields = {{90},{0,0.09:0.001:0.110};...
    {95},{0,0.09:0.001:0.110};...
    {100},{0,0.09:0.001:0.110}};

for pidx = 1:size(phifields,1);
    p1 = phifields{pidx,1};
    p1i = find(arrayfun(@(x)dblcmp(x,cell2mat(phifields{pidx,1})), phi1));
    
    % p2 = phifields{pidx,2};
    % fields_to_opt_struct = phifields{pidx,3};
    % p2i = find(arrayfun(@(x)dblcmp(x,cell2mat(phifields{pidx,2})), phi2));
    
    p2i = 1;
    fields_to_opt_struct = phifields{pidx,2};
    
    boandrdata{p1i,p2i} = cell(size(field,1),6);
    
    for fstrn = fields_to_opt_struct;
        k = find(arrayfun(@(x)dblcmp(x,fstrn), field));
        
        S = load(['s:\1Angle\',mols{imol,1},'-',num2str(fstrn),'VA.mat']);
        myexp = S.obj;

        thisdat = myexp.get_data(p1);
        % dat = myexp.get_data(p1,p2);
        
        dat = thisdat.indoOutput;
        dat = textscan(dat,'%s','delimiter','\n');
        dat = dat{1};
        dat = textscan(dat{2},'%s');
        dat = dat{1};

        efield = [str2double(dat{1}), str2double(dat{2}), str2double(dat{3})];

        thisdat.load_to_memory('indo','load');
        indo = thisdat.raw_indo;

        thisdat.generate_ampac_file('out', 'p:\');

        [rgs, rex, deltar, eszmat, ampac_energy, exindo, ~, numruns] = OptExcStateStructure('p:\', thisdat.ampac_hash,...
            'o:\', thisdat.indo_hash,...
            'field', efield,...
            'indo', indo,...
            'readifexist');

        % Second excited state
    %                 [rgs, rex, deltar, ~, numruns, eszmat, ampac_energy] = OptExcStateStructure('p:\', myexp.data(1).ampac_hash,...
    %                     'o:\', myexp.data(1).indo_hash,...
    %                     'field', efield,...
    %                     'indo', indo,...
    %                     'state', 3,...
    %                     'readifexist');

    %                 if (isempty(boandrdata))
    %                     natom = myexp.get_field('ampac.natom');
    %                     boandrdata = cell(natom,6);
    %                 end

        boandrdata{i,j}{k, 1} = fstrn;
        boandrdata{i,j}{k, 2} = rgs;
        boandrdata{i,j}{k, 3} = rex;
        boandrdata{i,j}{k, 4} = deltar;
        boandrdata{i,j}{k, 5} = numruns;
        boandrdata{i,j}{k, 6} = eszmat;
        % idx = idx + 1;

        exdeltaAM1{i,j}(k) = ampac_energy - myexp.data(1).Ehf;

        myexp.data(1).load_to_memory('ampac','load');
        amp = myexp.data(1).raw_ampac;
        amp.ampac_succeed = true;

        exindo = Indo.LoadExistingData(['o:\', myexp.data(1).indo_hash, '-new.ido'],[],[],[]);

        tmp = exindo.dipole(1,1);
        exdpgs{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

        tmp = exindo.dipole(2,2);
        exdp{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

        tmp = exindo.dipole(3,3);
        exdp2exc{i,j}(k) = sum(tmp .^ 2, 1) .^ (0.5);

        exnfill = exindo.nfilled;
        tmp = exindo.orbE([exnfill exnfill+1]);
        exhlgap{i,j}(k) = tmp(2) - tmp(1);

        tempds = ECEDataStruct(amp, exindo,[], [], [], [], []);
        exeexc{i,j}(k,1:25) = tempds.Eexc(:);

        exigs{i,j}(k) = exindo.esci(1);

        exopint{i,j}(k,1:25) = tempds.Tint(:);
        tmp = squeeze(exindo.r(2,:,:));
        for l = 1:25
            extpint{i,j}(k,l) = (exeexctemp(k,l) - exeexctemp(k,2)) * sum(tmp(l,:) .^ 2);
        end
    end
end

%% Save calcs for multiple molecules

save(['s:\',mols{imol,1},'-BOScript.mat']);

end
%% PLot GS Struct Energy vs. Field
figure(8)
hold on

plot_p1 = 175;
plot_p2 = [];

if (isempty(plot_p1))
    i = 1;
else
    i = find(arrayfun(@(x)dblcmp(x,plot_p1), phi1));
end

if (isempty(plot_p2))
    j = 1;
else
    j = find(arrayfun(@(x)dblcmp(x,plot_p2), phi2));
end


xaxis = field;

for k = 1:length(xaxis)
    plot(xaxis(k), eexc{i,j}(k,:) + igs{i,j}(k), 'g^');
end

maxopint = max(opint{i,j}(:));
for k = 1:length(xaxis)
    for l = 1:25
        plot(xaxis(k), eexc{i,j}(k,l) + igs{i,j}(k), 'bo', 'MarkerSize', (opint{i,j}(k,l)*30 / maxopint) + 1e-3);
    end
end

maxtpint = max(tpint{i,j}(:));

for k = 1:length(xaxis)
    for l = 1:25
        if (tpint{i,j}(k,l) > 0)
            plot(xaxis(k), eexc{i,j}(k,l) + igs{i,j}(k), 'rs', 'MarkerSize', (tpint{i,j}(k,l)*30 / maxtpint) + 1e-3);
        end
    end
end

%% Get slope and y-int of nBu. Calculate appx value of exciton binding energy
figure(8);  % Need to have energy vs. field plotted here already
[x,y] = ginput(3);
m = (y(2) - y(1)) / (x(2) - x(1));
b = y(1) - m*x(1);

disp(['The equation of the line is y = ', num2str(m), 'x - ', num2str(abs(b))]);

ebenergy = b - y(3);
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
    plot(xaxis(k), exeexc{i,j}(k,:) + igs{i,j}(k) + exdeltaAM1{i,j}(k), 'k^');
end

% maxexopint = max(exopint{i,j}(:));
for k = 1:length(xaxis)
    for l = 1:25
        plot(xaxis(k), exeexc{i,j}(k,l) + igs{i,j}(k) + exdeltaAM1{i,j}(k), 'mo', 'MarkerSize', (exopint{i,j}(k,l)*30 / maxopint) + 1e-3);
    end
end

% maxtpint = max(tpint{i,j}(:));

for k = 1:length(xaxis)
    for l = 1:25
        if (tpint{i,j}(k,l) > 0)
            plot(xaxis(k), exeexc{i,j}(k,l) + igs{i,j}(k) + exdeltaAM1{i,j}(k), 'cs', 'MarkerSize', (extpint{i,j}(k,l)*30 / maxtpint) + 1e-3);
        end
    end
end

% %% BEGIN RUNNING HERE
% end
%% Draw structure with delta r color coding

field_strength = 0.074;

S = load('s:\NoAngle\10-merPPV-0VA.mat');
myexp = S.obj;

myexp.data(1).generate_ampac_file('out');

[xyz, atom_type] = ampac_to_xyz(['p:\', myexp.data(1).ampac_hash, '.out']);
[exxyz, ~] = ampac_to_xyz(boandrdata{find(arrayfun(@(x)dblcmp(x, field_strength), field)), 6}); 
atomsbonds = get_connectivity(['p:\', myexp.data(1).ampac_hash, '.out']);

cidx = find(cellfun(@(x)strcmpi(x,'c'), atom_type));

% xyz = xyz(cidx,:);
% exxyz = exxyz(cidx,:);

figure(7);
hold on;
axis equal;
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2)) min(xyz(:,3))-1 max(xyz(:,3))+1]);
for i = 1:size(xyz,1)
    scatter3(xyz(i,1),xyz(i,2),xyz(i,3),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b');
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

plot_p1 = 90;
plot_p2 = [];

if (isempty(plot_p1))
    i = 1;
else
    i = find(arrayfun(@(x)dblcmp(x,plot_p1), phi1));
end

if (isempty(plot_p2))
    j = 1;
else
    j = find(arrayfun(@(x)dblcmp(x,plot_p2), phi2));
end

scatter(xaxis,reshape(dpgs{i,j}(1:length(xaxis)) .* 4.8,1,[]),'MarkerEdgeColor','g','Marker','*');
scatter(xaxis,reshape(dp{i,j}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','b','Marker','o');
scatter(xaxis,reshape(dp2exc{i,j}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','r','Marker','+');
scatter(xaxis,reshape(dp3exc{i,j}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','k','Marker','d');
scatter(xaxis,reshape(dp4exc{i,j}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','c','Marker','^');
% scatter(xaxis,reshape(dpmag{i,j}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','m','Marker','x');