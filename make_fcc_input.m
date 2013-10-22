%% Input files and parameters

% Files where the data are located
fn_gs_xyz = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-GSOpt-Modes.log';
fn_es_xyz = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-ESOpt.out';
fn_gs_norm_modes = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-GSOpt-Modes.log';
fn_es_norm_modes = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-ES-GetModes.out';
fn_gs_trdipole = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-GS-TD.log';
fn_es_trdipole = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene-ESOpt.out';

% File names to create for FCClasses input
gs_state_filename = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene.gs';
es_state_filename = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene.es';
trdp_out_filename = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene.td';
main_input_filename = 'C:\Users\clegaspi\Documents\Thiophenes\fcc\ethylene\ethylene.job';

% FCClasses input parameters
ex_st_of_interest = 1;  % This should probably always be 1, unless you care about emission from n=2 or something
calc_type = 'emi';      % Spectrum type. 'emi' = emission, 'abs' = absorption
ecd = 0;            % Do an ECD calculation? Yes=1, no=0
fc_ht = 'fc';       % Franck-Condon ('fc') or Hertzberg-Teller ('ht')

temperature = 100;  % Kelvin
boltz_pop_thresh = 0.3;     % minimum weight taken into account in the Boltzmann distribution.
                            % The weights is expressed as a fraction of the weight of the ground vibrational state (set to 1)
max_states_c1 = 30;     % Number of states to be considered for each oscillator in the computation of C1 class transitions
                        % (Overtones). Typical values 20-25
max_states_c2 = 25;     % Number of states to be computed for each oscillator in the computation of C2 class transitions
                        % (combinations of 2 modes). Typical values 15-20
acc_param = 1e6;        % Accuracy parameter. Nominal maximum number of integrals to be computed for each class
                        % Typical values : 1e6 - 1e7 (easy cases: rigid molecules with small displacements)
                        % 1e8 (standard)
                        % 1e9 - 1e10 (di±cult cases: very large molecules, and/or large displacements).
                        % Choices > 1e10 usually require very long computational times


% These parameters are for convolution. If we plan on doing our own convolution, ignore these                    
broadening_type = 'Gau';    % 'Gau' = Gaussian, 'Lor' = Lorentzian
en_plot_min = 1.4;      % Minimum energy for spectrum x-axis (eV)
en_plot_max = 2.4;      % Maximum energy for spectrum x-axis (eV)
npoints = 1001;         % Number of points in convolution
hwhm = 0.005;           % HWHM for lineshape (eV)



%% Read GS xyz

gsxyz = fileread(fn_gs_xyz);
gsxyz = textscan(gsxyz,'%s','delimiter','\n');
gsxyz = gsxyz{1};

found_opt = 0;
natoms = 0;
energy_idx = 0;

gs_energy = 0;
gs_geom = [];
gs_atom_type = [];

for i = 1:length(gsxyz)
    % Find the number of atoms if we haven't yet
    if (natoms == 0)
        if (~isempty(regexpi(gsxyz{i}, '.*NAtoms=\s*[0-9]+\s+.*')))
            dat = regexpi(gsxyz{i}, '.*NAtoms=\s*(?<natoms>[0-9]+)\s+.*', 'names');
            natoms = str2double(dat.natoms);
            gs_geom = zeros(natoms,3);
            gs_atom_type = zeros(natoms,1);
        end
    end
    
    % Have we gotten to the optimized section yet?
    if (~isempty(regexpi(gsxyz{i}, '.+Optimized Parameters.*')))
        found_opt = i;
    end
    
    % If in the optimized section, find the Cartesian coordinates in standard orientation
    % and read them in, along with the atomic number for creation of the
    % input file later.
    if (found_opt ~= 0)
        if (~isempty(regexpi(gsxyz{i}, '.*Standard orientation.*')))
            k = 1;
            for j = i+5:i+5+natoms-1
                tmp = regexpi(gsxyz{j}, '.*\d+\s+(?<atom_num>\d+)\s+\d+\s+(?<x>\S+)\s+(?<y>\S+)\s+(?<z>\S+)\s*', 'names');
                gs_atom_type(k) = str2double(tmp.atom_num);
                gs_geom(k,1) = str2double(tmp.x);
                gs_geom(k,2) = str2double(tmp.y);
                gs_geom(k,3) = str2double(tmp.z);
                k=k+1;
            end
        end
    end
    
    if (~isempty(regexpi(gsxyz{i}, '.*SCF Done:\s*E\(\S+\).*')))
        energy_idx = i;
    end
end

tmp = regexpi(gsxyz{energy_idx}, '.*SCF Done:\s*E\(\S+\) =(?<en>.+)A\.U\..*', 'names');
gs_energy = str2double(tmp.en);

%% Read GS norm modes

gsnorm = fileread(fn_gs_norm_modes);
gsnorm = textscan(gsnorm,'%s','delimiter','\n');
gsnorm = gsnorm{1};

if (~exist('natoms','var'))
    natoms = 0;
end
nmodes = 0;

gs_modes = [];
gs_freq = [];

for i = 1:length(gsnorm)    
    % Find the number of atoms if we haven't yet
    if (natoms == 0)
        if (~isempty(regexpi(gsnorm{i}, '.*NAtoms=\s*[0-9]+\s+.*')))
            dat = regexpi(gsnorm{i}, '.*NAtoms=\s*(?<natoms>[0-9]+)\s+.*', 'names');
            natoms = str2double(dat.natoms);
        end
    end
    
    % If we've found the normal modes...
    if (~isempty(regexpi(gsnorm{i}, '.*IR Intensities ---.*')))
        sidx = i-5;
        eidx = sidx;
        keep_going = true;
        
        % Figure out how many normal modes there are
        while (keep_going)
            if (isempty(regexpi(gsnorm{eidx+2}, '.*Frequencies ---.*')))
                keep_going = false;
                tmp = regexpi(gsnorm{eidx - 7 - 3*natoms}, '(?<nmodes>\d+)$', 'names');
                nmodes = str2double(tmp.nmodes);
                gs_modes = zeros(natoms,3,nmodes);
                gs_freq = zeros(nmodes,1);
            else
                eidx = eidx + 7 + 3*natoms;
            end
        end
        
        % Read them in as (atom number, xyz, mode number)
        for nk = 1:ceil(nmodes/5)
            for ni = 1:natoms
                for nj = 1:3
                    tmp = textscan(gsnorm{sidx+7+3*(ni-1)+(nj-1)},'%s');
                    tmp = tmp{1};
                    for idx = 4:length(tmp)
                        gs_modes(ni,nj,5*(nk-1)+idx-3) = str2double(tmp{idx});
                    end
                end
            end
            
            % Read in mode frequencies
            tmp = textscan(gsnorm{sidx+2},'%s');
            tmp = tmp{1};
            for idx = 3:length(tmp)
                gs_freq(5*(nk-1)+idx-2) = str2double(tmp{idx});
            end
            
            sidx = sidx + 7 + 3*natoms;
        end
    end
    
    % If we found everything, get out
    if (nmodes ~= 0)
        break
    end
end

%% GS Transition Dipoles

gstrdp = fileread(fn_gs_trdipole);
gstrdp = textscan(gstrdp,'%s','delimiter','\n');
gstrdp = gstrdp{1};
idx = 0;

gs_tr_dipole = zeros(3,1);

for i = 1:length(gstrdp)
    if (~isempty(regexpi(gstrdp{i}, '.*Ground to excited state transition electric dipole moments.*')))
        idx = i;
    end
end

idx = idx + 1 + ex_st_of_interest;
tmp = textscan(gstrdp{idx},'%s');
tmp = tmp{1};
for i = 1:3
    gs_tr_dipole(i) = str2double(tmp{i+1});
end

%% Read ES xyz

esxyz = fileread(fn_es_xyz);
esxyz = textscan(esxyz,'%s','delimiter','\n');
esxyz = esxyz{1};

found_opt = 0;
natoms = 0;
energy_idx = 0;

es_energy = 0;
es_geom = [];
es_atom_type = [];

for i = 1:length(esxyz)
    % Find the number of atoms if we haven't yet
    if (natoms == 0)
        if (~isempty(regexpi(esxyz{i}, '.*NAtoms=\s*[0-9]+\s+.*')))
            dat = regexpi(esxyz{i}, '.*NAtoms=\s*(?<natoms>[0-9]+)\s+.*', 'names');
            natoms = str2double(dat.natoms);
            es_geom = zeros(natoms,3);
            es_atom_type = zeros(natoms,1);
        end
    end
    
    % Have we gotten to the optimized section yet?
    if (~isempty(regexpi(esxyz{i}, '.+Optimized Parameters.*')))
        found_opt = i;
    end
    
    % If in the optimized section, find the Cartesian coordinates in standard orientation
    % and read them in, along with the atomic number for creation of the
    % input file later.
    if (found_opt ~= 0)
        if (~isempty(regexpi(esxyz{i}, '.*Standard orientation.*')))
            k = 1;
            for j = i+5:i+5+natoms-1
                tmp = regexpi(esxyz{j}, '.*\d+\s+(?<atom_num>\d+)\s+\d+\s+(?<x>\S+)\s+(?<y>\S+)\s+(?<z>\S+)\s*', 'names');
                es_atom_type(k) = str2double(tmp.atom_num);
                es_geom(k,1) = str2double(tmp.x);
                es_geom(k,2) = str2double(tmp.y);
                es_geom(k,3) = str2double(tmp.z);
                k=k+1;
            end
        end
    end
    
    if (~isempty(regexpi(esxyz{i}, '.*SCF Done:\s*E\(\S+\).*')))
        energy_idx = i;
    end
end

tmp = regexpi(esxyz{energy_idx}, '.*SCF Done:\s*E\(\S+\) =(?<en>.+)A\.U\..*', 'names');
es_energy = str2double(tmp.en);

%% Read ES norm modes

esnorm = fileread(fn_es_norm_modes);
esnorm = textscan(esnorm,'%s','delimiter','\n');
esnorm = esnorm{1};

if (~exist('natoms','var'))
    natoms = 0;
end
nmodes = 0;

es_modes = [];
es_freq = [];

for i = 1:length(esnorm)    
    % Find the number of atoms if we haven't yet
    if (natoms == 0)
        if (~isempty(regexpi(esnorm{i}, '.*NAtoms=\s*[0-9]+\s+.*')))
            dat = regexpi(esnorm{i}, '.*NAtoms=\s*(?<natoms>[0-9]+)\s+.*', 'names');
            natoms = str2double(dat.natoms);
        end
    end
    
    % If we've found the normal modes...
    if (~isempty(regexpi(esnorm{i}, '.*IR Intensities ---.*')))
        sidx = i-5;
        eidx = sidx;
        keep_going = true;
        
        % Figure out how many normal modes there are
        while (keep_going)
            if (isempty(regexpi(esnorm{eidx+2}, '.*Frequencies ---.*')))
                keep_going = false;
                tmp = regexpi(esnorm{eidx - 7 - 3*natoms}, '(?<nmodes>\d+)$', 'names');
                nmodes = str2double(tmp.nmodes);
                es_modes = zeros(natoms,3,nmodes);
                es_freq = zeros(nmodes,1);
            else
                eidx = eidx + 7 + 3*natoms;
            end
        end
        
        % Read them in as (atom number, xyz, mode number)
        for nk = 1:ceil(nmodes/5)
            for ni = 1:natoms
                for nj = 1:3
                    tmp = textscan(esnorm{sidx+7+3*(ni-1)+(nj-1)},'%s');
                    tmp = tmp{1};
                    for idx = 4:length(tmp)
                        es_modes(ni,nj,5*(nk-1)+idx-3) = str2double(tmp{idx});
                    end
                end
            end
            
            % Read in mode frequencies
            tmp = textscan(esnorm{sidx+2},'%s');
            tmp = tmp{1};
            for idx = 3:length(tmp)
                es_freq(5*(nk-1)+idx-2) = str2double(tmp{idx});
            end
            
            sidx = sidx + 7 + 3*natoms;
        end
    end
    
    % If we found everything, get out
    if (nmodes ~= 0)
        break
    end
end

%% ES Transition Dipoles

estrdp = fileread(fn_es_trdipole);
estrdp = textscan(estrdp,'%s','delimiter','\n');
estrdp = estrdp{1};
idx = 0;

es_tr_dipole = zeros(3,1);

for i = 1:length(estrdp)
    if (~isempty(regexpi(estrdp{i}, '.*Ground to excited state transition electric dipole moments.*')))
        idx = i;
    end
end

idx = idx + 1 + ex_st_of_interest;
tmp = textscan(estrdp{idx},'%s');
tmp = tmp{1};
for i = 1:3
    es_tr_dipole(i) = str2double(tmp{i+1});
end


%% Output GS FCClasses File (F01 or F02)

% This is the "state" file for the ground state
gs_fid = fopen(gs_state_filename,'w');

% Cartesian coordinates in standard orientation
for i = 1:natoms
    for j = 1:3
        fprintf(gs_fid,'%E\n',gs_geom(i,j));
    end
end

% Normal modes
for i = 1:natoms
    for j = 1:3
        for k = 1:nmodes
            fprintf(gs_fid, '%E\n', gs_modes(i,j,k));
        end
    end
end

% Frequencies of normal modes
for i = 1:nmodes
    fprintf(gs_fid, '%E\n', gs_freq(i));
end
            
fclose(gs_fid);

%% Output ES FCClasses File (F01 or F02)

% This is the "state" file for the excited state
es_fid = fopen(es_state_filename, 'w');

% Cartesian coordinates in standard orientation
for i = 1:natoms
    for j = 1:3
        fprintf(es_fid,'%E\n',es_geom(i,j));
    end
end

% Normal modes
for i = 1:natoms
    for j = 1:3
        for k = 1:nmodes
            fprintf(es_fid, '%E\n', es_modes(i,j,k));
        end
    end
end

% Frequencies of normal modes
for i = 1:nmodes
    fprintf(es_fid, '%E\n', es_freq(i));
end
            
fclose(es_fid);

%% Output Transition dipole file (F03)

fid = fopen(trdp_out_filename, 'w');

fprintf(fid, '%E\t%E\t%E\n', gs_tr_dipole(1:3));
fprintf(fid, '%E\t%E\t%E', es_tr_dipole(1:3));

fclose(fid);

%% Output main input file for FCClasses

fid = fopen(main_input_filename, 'w');

fprintf(fid, '%u\t\t\t Number of atoms\n', natoms);
fprintf(fid, '%u\t\t\t Number of normal modes\n', nmodes);

for i = 1:natoms
    mass = get_atomic_mass(gs_atom_type(i));   % ES or GS, doesn't matter
    fprintf(fid, '%3.4f\t\t\t Mass of atom %u\n', mass, i);
end

% en_diff = abs(es_energy - gs_energy) * 27.2107;     % Hartrees to eV
en_diff = 3.3132;   % 2T
fprintf(fid, '%E\t\t\t Adiabatic energy difference\n', en_diff);

if (ecd == 0)
    tmp = 'ECDNO';
else
    tmp = 'ECDYES';
end
fprintf(fid, '''%s'' ''%s'' ''%s''\t\t\t Calculation type\n', calc_type, tmp, fc_ht);

fprintf(fid, '%E %E\t\t\t Temperature and Boltzmann\n''D''\t\t\t Option for temp calculation use, use D\n',...
    temperature, boltz_pop_thresh);

if (strcmpi(calc_type,'emi'))
    fprintf(fid, '''%s''\t\t\t Originating state info\n', es_state_filename);
    fprintf(fid, '''%s''\t\t\t Ending state info\n', gs_state_filename);
else
    fprintf(fid, '''%s''\t\t\t Originating state info\n', gs_state_filename);
    fprintf(fid, '''%s''\t\t\t Ending state info\n', es_state_filename);
end

fprintf(fid, '''%s''\t\t\t Transition dipole info\n', trdp_out_filename);
fprintf(fid, '''%s''\t\t\t Magnetic dipole info\n', 'F04 file goes here if we plan to use eventually');
fprintf(fid, '0\t\t\t Option for rotation\n');
fprintf(fid, '%u\t\t\t Max states in C1\n', max_states_c1);
fprintf(fid, '%u\t\t\t Max states in C2\n', max_states_c2);
fprintf(fid, '%E\t\t\t Accuracy parameter\n', acc_param);
fprintf(fid, '''%s'' %E %E %u %E\t\t\t Convolution parameters\n',...
    broadening_type, en_plot_min, en_plot_max, npoints, hwhm);

fclose(fid);
            