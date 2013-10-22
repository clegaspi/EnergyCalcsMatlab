%% Paths and settings

cygbash = 'C:\cygwin\bin\bash.exe';
homedir = 'C:\cygwin\home\clegaspi';
username = 'clegaspi';
server = 'c60.chem.cmu.edu';
sqmpath = '~/amber12/AmberTools/bin/sqm';    % Case-sensitive
jobname = 'bithiophene';
%% Running a loop for testing
% [xyz, atomtype] = ampac_to_xyz('C:\Users\Christian\Documents\Research\Yaron\dyes2\data\8-merPPV\Exp\NoAngle\8-merPPV.out');
[xyz, atomtype] = ampac_to_xyz('C:\Users\clegaspi\Documents\MATLAB\data\trans-bithiophene\trans-bithiophene-slighttwist.out');
efv = xyz(4,:)-xyz(6,:);
efv = efv / norm(efv);

for field = 0
    jobname = ['trans-bithiophene-',num2str(field),'VA'];

%% Write job file
fid = fopen([homedir,'\',jobname,'.in'],'w');

fprintf(fid, '%s\n', jobname);
fprintf(fid, '%s\n', ' &qmmm');
fprintf(fid, '\t%s\n', 'qm_theory = ''AM1'',');
fprintf(fid, '\t%s\n', 'exst_method = 1,');
fprintf(fid, '\t%s\n', 'qmqmdx = 2,');
fprintf(fid, '\t%s\n', 'qmcharge = 0,');
fprintf(fid, '\t%s\n', 'scfconv = 1.0d-10,');
fprintf(fid, '\t%s\n', 'errconv = 1.0d-6,');
fprintf(fid, '\t%s\n', 'tight_p_conv = 0,');
fprintf(fid, '\t%s\n', 'diag_routine = 2,');
fprintf(fid, '\t%s\n', 'fock_predict = 1,');
fprintf(fid, '\t%s\n', 'itrmax = 100,');
fprintf(fid, '\t%s\n', 'efield_flag = 1,');
fprintf(fid, '\t%s\n', ['efield_vector = ',num2str(efv(1)*field,'%5.9f'),' ',num2str(efv(2)*field,'%5.9f'),...
    ' ',num2str(efv(3)*field,'%5.9f'),',']);
fprintf(fid, '\t%s\n', 'maxcyc = 0,');
fprintf(fid, '\t%s\n', 'grms_tol = 0.02,');
fprintf(fid, '\t%s\n', 'ntpr = 101,');
fprintf(fid, '\t%s\n', 'verbosity = 5,');
fprintf(fid, '\t%s\n', 'excN = 25,');
fprintf(fid, '\t%s\n', 'struct_opt_state = 1,');
fprintf(fid, '%s\n', ' /');

for i = 1:length(atomtype)
    if (strcmpi(atomtype{i},'H'))
        atomic_num = 1;
    elseif (strcmpi(atomtype{i},'C'))
        atomic_num = 6;
    elseif (strcmpi(atomtype{i},'O'))
        atomic_num = 8;
    elseif (strcmpi(atomtype{i},'N'))
        atomic_num = 7;
    elseif (strcmpi(atomtype{i},'S'))
        atomic_num = 16;
    end
    fprintf(fid, ' %i %s\t%5.6f\t%5.6f\t%5.6f\n', atomic_num, atomtype{i}, xyz(i,1), xyz(i,2), xyz(i,3));
end

fclose(fid);

%% Write local bash script
fid = fopen([homedir,'\local_script.sh'],'w');

fprintf(fid, '%s\n', '#/bin/bash');
fprintf(fid, '%s\n', ['sftp -b - ', username, '@', server, ' <<EOF']);
fprintf(fid, '%s\n', ['put ',jobname,'.in sqm_jobs/']);
fprintf(fid, '%s\n', 'bye');
fprintf(fid, '%s\n', 'EOF');
fprintf(fid, '%s\n', ['ssh ', username, '@', server, ' ''/bin/bash -s'' < ~/remote_script.sh']);
fprintf(fid, '%s\n', ['sftp -b - ', username, '@', server, ' <<EOF']);
fprintf(fid, '%s\n', ['get sqm_jobs/',jobname,'.out .']);
fprintf(fid, '%s\n', 'bye');
fprintf(fid, '%s\n', 'EOF');

fclose(fid);

%% Write remote bash script
fid = fopen([homedir,'\remote_script.sh'],'w');

fprintf(fid, '%s\n', '#/bin/bash');
fprintf(fid, '%s\n', 'cd ~/sqm_jobs');
fprintf(fid, '%s\n', [sqmpath,' -i ',jobname,'.in -o ',jobname,'.out -O']);

fclose(fid);

%% Execute local bash script

system([cygbash, ' -li ~/local_script.sh']);

%% Terminate test loop

end

%% Interpret results

field = 0:0.001:0.150;
osc = zeros(length(field),25);
eexc = zeros(length(field),25);
gsen = zeros(1,length(field));
dpgs = zeros(1,length(field));
dp = zeros(length(field),25);

for i = 1:length(field)
    if (~exist([homedir,'\8-merPPV-',num2str(field(i)),'VA.out'],'file'))
        disp(['Did not find field ',num2str(field(i))]);
        keyboard;
    end
    mytext = fileread([homedir,'\8-merPPV-',num2str(field(i)),'VA.out']);
    mytext = textscan(mytext, '%s', 'delimiter', '\n');
    mytext = mytext{1};
    
    idx = find(cellfun(@(x)strcmp(x,'Total energy of ground state (eV)'),mytext))+1;
    disp(i);
    gsen(i) = str2double(mytext{idx});
    
    idx = find(cellfun(@(x)strcmp(x,'Frequencies (eV) and Oscillator strengths (unitless)'),mytext))+1;
    for j = 1:25
        tmp = textscan(mytext{idx+j},'%s');
        tmp = tmp{1};
        eexc(i,j) = str2double(tmp{2});
        osc(i,j) = str2double(tmp{6});
    end
    
    idx = find(cellfun(@(x)strcmp(x,'DIPOLES OF EXCITED STATES'),mytext))+3;
    tmp = textscan(mytext{idx-5},'%s');
    tmp = tmp{1};
    dpgs(i) = str2double(tmp{6});
    for j = 1:25
        tmp = textscan(mytext{idx+(j-1)*3},'%s');
        tmp = tmp{1};
        dp(i,j) = str2double(tmp{6});
    end
end

%% GS Structure Energy vs. Field
figure(8)
hold on

xaxis = field;

for k = 1:length(xaxis)
    plot(xaxis(k), gsen(k), 'g^');
    plot(xaxis(k), eexc(k,:) + gsen(k), 'g^');
    % plot(xaxis(k), eexc(k,:), 'g^');
end

maxopint = max(osc(:));
for k = 1:length(xaxis)
    for l = 1:25
        plot(xaxis(k), eexc(k,l) + gsen(k), 'bo', 'MarkerSize', (osc(k,l)*30 / maxopint) + 1e-3);
        % plot(xaxis(k), eexc(k,l), 'bo', 'MarkerSize', (osc(k,l)*30 / maxopint) + 1e-3);
    end
end

% maxtpint = max(tpint{1,1}(:));

% for k = 1:length(xaxis)
%     for l = 1:25
%         if (tpint{1,1}(k,l) > 0)
%             % plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k), 'rs', 'MarkerSize', (tpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
%             plot(xaxis(k), eexc{1,1}(k,l), 'rs', 'MarkerSize', (tpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
%         end
%     end
% end

%% Get slope and y-int of nBu. Calculate appx value of exciton binding energy
figure(8);  % Need to have energy vs. field plotted here already
[x,y] = ginput(3); % Click two points on the sloped line that crosses 1Bu and once on the 1Bu
m = (y(2) - y(1)) / (x(2) - x(1));
b = y(1) - m*x(1);

disp(['The equation of the line is y = ', num2str(m), 'x - ', num2str(abs(b))]);

ebenergy = b - y(3);
disp(['The exciton binding energy is approximately ', num2str(ebenergy), ' eV']);

%% Dipole moments
figure(17);
hold on;
xaxis = field;

% Dipoles from SQM are output in Debye
scatter(xaxis,dpgs(1:length(xaxis)),'MarkerEdgeColor','g','Marker','*');
scatter(xaxis,reshape(dp(1:length(xaxis),1),1,[]),'MarkerEdgeColor','b','Marker','o');
scatter(xaxis,reshape(dp(1:length(xaxis),2),1,[]),'MarkerEdgeColor','r','Marker','+');
% scatter(xaxis,reshape(dp3exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','k','Marker','d');
% scatter(xaxis,reshape(dp4exc{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','c','Marker','^');
% scatter(xaxis,reshape(dpmag{k}(1:length(xaxis)).* 4.8,1,[]),'MarkerEdgeColor','m','Marker','x');

%% Delta R Color Coding

%% File params
testdir = homedir;
dir = '8-merPPV';
df = 'geomopt-0.140VA.out';
mol2file = '8-merPPV.mol2';

%% Read coordinates
datfile = fileread([testdir,'\', df]);
datfile = textscan(datfile,'%s','delimiter','\n');
datfile = datfile{1};

idx = find(cellfun(@(x)strcmp(x,'QMMM: QM Region Cartesian Coordinates (*=link atom) '), datfile))+2;

% natoms = length(init);
% gs = zeros(natoms,3);
% es = gs;
i = 1;
while (~isempty(datfile{idx(1)}))
    tmp = textscan(datfile{idx(1)},'%s');
    tmp = tmp{1};
    gs(i,1:3) = cellfun(@(x)str2double(x), tmp(5:7))';
    
    tmp = textscan(datfile{idx(2)},'%s');
    tmp = tmp{1};
    es(i,1:3) = cellfun(@(x)str2double(x), tmp(5:7))';
    
    idx = idx + 1;
    i = i + 1;
end

natoms = size(gs,1);

%% read bond structure and atom types

mf = fileread([testdir,'\',mol2file]);

mf = textscan(mf,'%s','delimiter','\n');
mf = mf{1};

bidx = find(cellfun(@(x)strcmp(x, '@<TRIPOS>BOND'), mf)) + 1;
zidx = find(cellfun(@(x)strcmp(x, '@<TRIPOS>ATOM'), mf)) + 1;

zmf = {mf{zidx:bidx-1}};
bmf = {mf{bidx:end}};

bmf = cellfun(@(x)textscan(x,'%s'), bmf);
zmf = cellfun(@(x)textscan(x,'%s'), zmf);

atoms = zeros(1,natoms);
for i = 1:natoms
    switch zmf{i}{end}
        case 'C'
            atoms(i) = 6;
        case 'H'
            atoms(i) = 1;
        case 'O'
            atoms(i) = 8;
        case 'N'
            atoms(i) = 7;
        case 'F'
            atoms(i) = 9;
    end
end

bonds = zeros(length(bmf),2);

for i = 1:length(bmf)
    bonds(i,:) = cellfun(@(x)str2double(x), bmf{i}(2:3));
end

%% Draw atoms

figure(1)
hold on
axis([min(gs(:,1))-1 max(gs(:,1))+1 min(gs(:,2))-1 max(gs(:,2))+1 min(gs(:,3)) max(gs(:,3))]);
axis equal;
whitebg([119 136 153] / 255);
for i=1:natoms
    if (atoms(i) == 6)
        color = 'k';
    elseif (atoms(i) == 1)
        color = 'w';
    elseif (atoms(i) == 8)
        color = 'r';
    elseif (atoms(i) == 7)
        color = 'b';
    elseif (atoms(i) == 9)
        color = 'y';
    end
    scatter3(gs(i,1),gs(i,2),gs(i,3),[color,'o'],'LineWidth',5);
end

%% Calculate delta x and draw bonds onto GS structure
deltax = zeros(1,length(bonds));
for i = 1:length(bonds)
    deltax(i) = sqrt(sum((es(bonds(i,1),:) - es(bonds(i,2),:)).^2))...
        - sqrt(sum((gs(bonds(i,1),:) - gs(bonds(i,2),:)).^2));
end
for i = 1:length(bonds)
    if (deltax(i) < 0)
        color = 'r';
    elseif (deltax(i) > 0)
        color = 'c';
    else
        color = 'k';
    end
    dat = [gs(bonds(i,1),:)' gs(bonds(i,2),:)'];
    line(dat(1,:), dat(2,:), dat(3,:),'color',color,'linewidth',(abs(deltax(i)) * 5 / max(abs(deltax))) + 1e-5);
end

%% Bond Length Alternation Plot
% Enter atom numbers of S-D-S units where 1-2 and 3-4 are singles and 2-3
% is a double
bla_anums = [3 7 8 9; 14:17; 22:25; 30:33; 38:41; 46:49; 54:57];

gs_bla = zeros(size(bla_anums, 1),1);
ex_bla = gs_bla;

for i = 1:length(gs_bla)
    sb1 = sum((gs(bla_anums(i,1),:) - gs(bla_anums(i,2),:)) .^ 2) ^ 0.5;
    sb2 = sum((gs(bla_anums(i,3),:) - gs(bla_anums(i,4),:)) .^ 2) ^ 0.5;
    db = sum((gs(bla_anums(i,2),:) - gs(bla_anums(i,3),:)) .^ 2) ^ 0.5;
    gs_bla(i) = ((sb1 + sb2) / 2) - db;
    
    sb1 = sum((es(bla_anums(i,1),:) - es(bla_anums(i,2),:)) .^ 2) ^ 0.5;
    sb2 = sum((es(bla_anums(i,3),:) - es(bla_anums(i,4),:)) .^ 2) ^ 0.5;
    db = sum((es(bla_anums(i,2),:) - es(bla_anums(i,3),:)) .^ 2) ^ 0.5;
    ex_bla(i) = ((sb1 + sb2) / 2) - db;
end

figure(5);
hold on
plot(1:length(gs_bla), gs_bla, '-g');
plot(1:length(ex_bla), ex_bla, '-r');