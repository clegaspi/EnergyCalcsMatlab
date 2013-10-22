%% Read in optimized ground state z-matrix
function [rgs, rex, deltar, newzmatrix, ampac_energy, indo, varargout] = ...
    OptExcStateStructure(ampac_pathonly, ampac_nameonly, indo_pathonly, indo_nameonly, varargin)
    
% ampac_pathonly = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\DMG-8mer\';
% ampac_nameonly = '8-merPPVampac';
    
    % EXECUTABLE PATHWAYS, NO QUOTES
    AmpacEXE = 'C:\Program Files\Semichem, Inc\Ampac-10.1\ampac.exe';
    OpenBabelEXE = 'C:\Program Files (x86)\OpenBabel-2.3.2\babel.exe';
    
    
    indoidx = find(cellfun(@(x)isequal(lower(x),'indo'), varargin));
    efieldidx = find(cellfun(@(x)isequal(lower(x),'field'), varargin));
    stateidx = find(cellfun(@(x)isequal(lower(x),'state'), varargin));
    read_if_exist = any(cellfun(@(x)isequal(lower(x),'readifexist'), varargin));
    blalgidx = find(cellfun(@(x)isequal(lower(x),'algorithm'), varargin));
    ppidx = find(cellfun(@(x)isequal(lower(x),'c'), varargin));
    nstates = find(cellfun(@(x)isequal(lower(x),'nstates'), varargin));     % How many states to calculate in INDO
    
    if (~isempty(blalgidx))
        blbo_algorithm = varargin{blalgidx+1};
    else
        blbo_algorithm = 'pauling';
    end
    
    if (~isempty(ppidx))
        pauling_param = varargin{ppidx+1};
    else
        pauling_param = 0.3;
    end
    
    if (~isempty(stateidx))
        state_to_opt = varargin{stateidx+1};
    else
        state_to_opt = 2;   % First excited state
    end
    
    if (~isempty(nstates))
        nstates = varargin{nstates+1};
    else
        nstates = state_to_opt + 1; % The +1 is just a formality...
    end
    
    ampac_filepath = [ampac_pathonly, ampac_nameonly];
    indo_filepath = [indo_pathonly, indo_nameonly];
    
    exist_dont_overwrite = false;
    
    if (read_if_exist)
        if (exist([indo_filepath, '-new.ido'], 'file') && exist([indo_filepath, '-new-dm.bin'], 'file'))
            exist_dont_overwrite = true;
        end
    end
    
    copyfile([ampac_filepath, '.out'],[ampac_filepath, '-new.out']);
    
    if (~exist_dont_overwrite)
        
        copyfile([indo_filepath, '-dm.bin'],[indo_filepath, '-new-dm.bin']);
    end
        
    if (~isempty(indoidx))
        if (exist_dont_overwrite)
            indo = Indo.LoadExistingData([indo_filepath,'-new.ido'],[],[],[]);
        else
            indo = varargin{indoidx+1};
        end
        indo_nameonly = [indo_nameonly, '-new'];
        indo_filepath = [indo_pathonly, indo_nameonly];
    else
        if (~exist_dont_overwrite)
            copyfile([indo_filepath, '.ido'],[indo_filepath, '-new.ido']);
        end
        indo_nameonly = [indo_nameonly, '-new'];
        indo_filepath = [indo_pathonly, indo_nameonly];
        indo = Indo.LoadExistingData([indo_filepath,'.ido'],[],[],[]);
    end
        
       
    ampac_nameonly = [ampac_nameonly, '-new'];
    ampac_filepath = [ampac_pathonly, ampac_nameonly];
    
    
    getOut = false;
    oldbls = [];
    numruns = 0;
    nrunsincreasing = 0;
    lowest_diff = [];

%    while (~getOut)
    
    zmatrix = ampac_to_zmatrix([ampac_filepath, '.out']);
    
    newzmatrix = zmatrix;
    
    [atomsbonds, mf] = get_connectivity(zmatrix);
    
    natoms = size(zmatrix, 1);

    %% Get rings

    rings = {};     % known rings
    ringset = [];   % current path we're following
    atomids = [0];  % atom numbers of known ring atoms
    gotaring = false;   % bool to let us get out of our current path if we've found a ring

    for i = 1:length(atomsbonds)    % loop over all atom numbers
        if (~ismember(i, atomids))  % check to see if atom is one of the known ring atoms
            ringset(end+1) = i;     % First atom in path
            for j = 1:length(atomsbonds{i}) % Loop over all paths from atom i
                if (length(atomsbonds{atomsbonds{i}(j)}) > 1)   % Check for backtracking
                    ringset(end+1) = atomsbonds{i}(j);  % Select a path from i
                    if (~gotaring)  % Keep going if we're not trying to get out from a found ring
                        for k = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom j
                            if (atomsbonds{ringset(end)}(k) ~= i && ~gotaring)  % Check for backtracking & getting out
                                ringset(end+1) = atomsbonds{ringset(end)}(k); % Select path from j
                                for l = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom k
                                    if (atomsbonds{ringset(end)}(l) ~= ringset(end-1) && ~gotaring) % Check for backtracking
                                        ringset(end+1) = atomsbonds{ringset(end)}(l);   % Select path from k
                                        for m = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom l
                                            if (atomsbonds{ringset(end)}(m) ~= ringset(end-1) && ~gotaring) % Check for backtracking
                                                ringset(end+1) = atomsbonds{ringset(end)}(m);   % Select path from l
                                                for n = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom m
                                                    if (atomsbonds{ringset(end)}(n) ~= ringset(end-1))  % Check for backtracking
                                                        if (~ismember(i, atomsbonds{ringset(end)})) % Check to see if any path fom m is to atom i
                                                            if (~gotaring)  % If not, check to get out and continue
                                                                    ringset(end+1) = atomsbonds{ringset(end)}(n); % Select a path from m
                                                                    if (ismember(i, atomsbonds{ringset(end)}) && atomsbonds{ringset(end)}(n) ~= ringset(end-1)) % Check if any path from atom n is to atom i
                                                                        rings{end+1} = ringset; % If so, we found a ring, record it
                                                                        atomids = unique([atomids ringset]);
                                                                        gotaring = true; % and get out
                                                                    end
                                                                    ringset = ringset(1:end-1); % remove n from current path
                                                            end
                                                        elseif (~gotaring)  % if we're not trying to get out and we found a 5-membered ring, record it and get out
                                                            rings{end+1} = ringset;
                                                            atomids = unique([atomids ringset]);
                                                            gotaring = true;
                                                        end
                                                    end
                                                end 
                                                ringset = ringset(1:end-1); % remove m from current path
                                            end
                                        end
                                        ringset = ringset(1:end-1); % remove l from current path
                                    end
                                end
                                ringset = ringset(1:end-1); % remove k from current path
                            end
                        end

                    end
                    ringset = ringset(1:end-1); % remove j from current path
                end
            end
        end
        ringset = ringset(1:end-1); % remove i from current path
        gotaring = false; % start looking again, excluding atoms we marked as in a ring already
    end

    %% Get z-matrix angles

    ringonly = [];
    ringinc = [];
    noring = [];

    for i = 1:natoms-2
        aset = [str2double(zmatrix{i+2,1}) str2double(zmatrix{i+2,6}) str2double(zmatrix{i+2,7})];
        switch (length(intersect(aset, atomids)))
            case 3
                ringonly(end+1,:) = aset;
            case 2
                ringinc(end+1,:) = aset;
            otherwise
                noring(end+1,:) = aset;          
        end
    end
    
    ring_info = struct('rings',rings,'ringonly',ringonly,'ringinc',ringinc);

    % allbanums = [ringonly; ringinc; noring];

    %% Get all C-C and C-H bonds

    CCbonds = {};
    CHbonds = {};

    for i = 1:length(mf)
        if (strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'c') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'c'))
            myset = sort([str2double(mf{i}(2)) str2double(mf{i}(3))]);
            if (~any(cellfun(@(x)isequal(x, myset), CCbonds)))
                CCbonds{end+1} = myset;
            end
        elseif ((strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'c') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'h')) || ...
                (strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'h') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'c')))
            myset = sort([str2double(mf{i}(2)) str2double(mf{i}(3))]);
            if (~any(cellfun(@(x)isequal(x, myset), CHbonds)))
                CHbonds{end+1} = myset;
            end
        end
    end

    temp = reshape(cell2mat(CCbonds),2, length(CCbonds))';
    [b,ix] = sort(temp,1);
    CCbonds = [b(:,1), temp(ix(:,1), 2)];

    temp = reshape(cell2mat(CHbonds),2, length(CHbonds))';
    [b,ix] = sort(temp,1);
    CHbonds = [b(:,1), temp(ix(:,1), 2)];
    
    
    oldbls = cellfun(@(x)str2double(x),{zmatrix{2:end,3}});
    rgs = oldbls;
    rex = rgs;
    gsbo = [];
    gsbobyatom = [];

    while (~getOut)
        numruns = numruns + 1;
        % newzmatrix = zmatrix;
        
        %% Read in ground state density matrix file

        fid = fopen([indo_filepath,'-dm.bin']);

        matdim = fread(fid,1,'int');

        dm = fread(fid,[matdim,matdim],'double');

        fclose(fid);
        
        %% Read GS density matrix and calculate GS C-C bond orders
        
        if (isempty(gsbo) && ~exist_dont_overwrite)
            gsbo = zeros(size(CCbonds,1),1);
            gsbobyatom = zeros(natoms);

            for i = 1:size(CCbonds,1)
                gsblmat = dm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
                gsbo(i) = sum(gsblmat(:));
                gsbobyatom(CCbonds(i,1), CCbonds(i,2)) = gsbo(i);
                gsbobyatom(CCbonds(i,2), CCbonds(i,1)) = gsbo(i);
            end
        elseif (exist_dont_overwrite)
            fid = fopen([indo_filepath(1:end-4),'-dm.bin']);
            matdim = fread(fid,1,'int');
            gsdm = fread(fid,[matdim,matdim],'double');
            fclose(fid);
            
            gsbo = zeros(size(CCbonds,1),1);
            gsbobyatom = zeros(natoms);
            
            for i = 1:size(CCbonds,1)
                gsblmat = gsdm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
                gsbo(i) = sum(gsblmat(:));
                gsbobyatom(CCbonds(i,1), CCbonds(i,2)) = gsbo(i);
                gsbobyatom(CCbonds(i,2), CCbonds(i,1)) = gsbo(i);
            end
        end
        
        %% Calculate DM change and generate excited state DM
        
        deltadm = diffDensity(indo, state_to_opt);
        
        dm = dm + deltadm;

        %% Calculate excited state bond order of C-C bonds

        bo = zeros(size(CCbonds,1),1);
        bobyatom = zeros(natoms);

        for i = 1:size(CCbonds,1)
            blmat = dm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
            bo(i) = sum(blmat(:));
            bobyatom(CCbonds(i,1), CCbonds(i,2)) = bo(i);
            bobyatom(CCbonds(i,2), CCbonds(i,1)) = bo(i);
        end
        
        % Columns 1-3: Atom 1, Atom 2, Calculated ES Bond Order
        % atomsandbo = [CCbonds(:,1) CCbonds(:,2) bo];

        %% Calculate new bond lengths and ring angles

        
        % Apply Pauling bond order equation to generate delta(r) for new
        % bond lengths
        
        newblsbyatom = zeros(natoms,natoms);
        [xyz, ~] = ampac_to_xyz(zmatrix);
        
        for i = 1:size(atomsbonds,1)
            for j = atomsbonds{i};
                trgs = gsbobyatom(i,j);
                trex = bobyatom(i,j);
                deltar = pauling_param * log(trgs / trex);
                
%                 if (numruns < 5)
%                     deltar = deltar * 2;
%                 end

                oldbl = sum((xyz(i,:) - xyz(j,:)) .^ 2) ^ 0.5;
                newbl = oldbl + deltar;

                newblsbyatom(i,j) = newbl;
                newblsbyatom(j,i) = newbl;
            end
        end
        
        available_algorithms = {'pauling','paulingwitham1'};
        switch (find(cellfun(@(x)strcmpi(blbo_algorithm, x), available_algorithms)))
            case 1
                newzmatrix = distort_geometry(zmatrix, ring_info, newblsbyatom);
            case 2
                newzmatrix = distort_geometry(zmatrix, ring_info, newblsbyatom, 'am1optimize');
            otherwise
                throw(MException('OptExcStateStructure:bad_algorithm','The BLBO algorithm specified does not exist!'));
        end
            
        
        
        %% Get current bond lengths and check them
        
        temp = cellfun(@(x)str2double(x),{newzmatrix{2:end,3}});
        diff = max(abs(temp - rex));
        rex = temp;
        
        if (~exist_dont_overwrite)
            disp(['OptExcStStruct: run ', num2str(numruns), ', diff ', num2str(diff)]);
        end
        
        if (diff < 1e-4 || exist_dont_overwrite || nrunsincreasing > 9)
            if (nrunsincreasing > 9)
                bobyatom = save_bobyatom;
                rex = save_rex;
                newzmatrix = save_newzmatrix;
                disp(['Structure optimization halted at diff = ',num2str(lowest_diff),'. Outputting those results.']);
            end
            
            getOut = true;
            % deltar = rex - rgs;
            deltar = zeros(natoms-1,1);
            atom_nums_in_bonds = zeros(natoms-1,2);
            
            for i = 1:natoms-1
                trgs = gsbobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                trex = bobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                atom_nums_in_bonds(i,1:2) = [str2double(zmatrix{i+1,1}), str2double(zmatrix{i+1,6})];
                deltar(i) = pauling_param * log(trgs / trex);
            end
            
            rgs = rgs';
            rex = rex';
            
            if (nargout > 6)
                varargout{1} = ring_info;
            end
            if (nargout > 7)
                varargout{2} = numruns;
            end
            if (nargout > 8)
                varargout{3} = atom_nums_in_bonds;
            end
        end
        
        if (isempty(lowest_diff) || diff < lowest_diff)
            lowest_diff = diff;
            save_bobyatom = bobyatom;
            save_rex = rex;
            save_newzmatrix = newzmatrix;
            nrunsincreasing = 0;
        else
            nrunsincreasing = nrunsincreasing + 1;
        end
            
        
        %% Write out new z-matrix as Ampac DAT for single-point calculation

        zmatrix_to_ampac(newzmatrix, ampac_pathonly, ampac_nameonly);


        %% Run single point calculation to generate OUT file and then get Ampac energy

        [~,~] = system(['"', AmpacEXE, '" "', ampac_filepath, '.dat"']);
        
        es_ampac = parseAmpac(ampac_filepath);
        
        ampac_energy = es_ampac.Hf / 23.05; % Convert to eV
        
        %% Send new structure to INDO if we need to
        
        if (~getOut)            
            if (~isempty(efieldidx))
                efield = varargin{efieldidx+1};
            else
                efield = [0 0 0];
            end
            
            res = [];
            res.charge = 0;
            res.norbs = 500;
            res.nstates = nstates;
            res.field = efield;
            res.initial_shiftc = 82.0;
            res.initial_shift_step = 1.0;
            res.min_shift_step = 0.1;
            res.max_shift_step = 10.0;
            res.initial_second_shift_step = 0.5;
            res.min_second_shift_step = 0.01;
            res.max_second_shift_step = 0.5;
            res.initial_eeint = 1.0;
            res.initial_eestep = 0.0;
            res.min_eestep = 0.0;
            res.max_eestep = 0.0;
            res.initial_conv = 1e-3;
            res.min_conv = 1e-10;
            res.max_inner_iter = 5000;
            res.max_iter = 150000;
            res.dm_guess = [indo_filepath, '-dm.bin'];
            res.try_default_first = true;
            res.output_dm = true;
            res.pot_file = [];

            indo = Indo(res, ampac_pathonly, ampac_nameonly);
            if (~strcmp(ampac_filepath, indo_filepath))
                movefile([ampac_pathonly, ampac_nameonly, '.ido'], [indo_pathonly, indo_nameonly, '.ido']);
                movefile([ampac_pathonly, ampac_nameonly, '-dm.bin'], [indo_pathonly, indo_nameonly, '-dm.bin']);
            end
            % delete([ampac_pathonly, ampac_nameonly, '.ipf']);
        end
    end
end