classdef Indo < handle
    properties (SetAccess = private)
        % Input parameters
        config     % see defaultConfig() for contents
        dataPathIn   % directory for the data (do not end with \)
        jobNameIn    % will read structure from ampac out file  jobname.out
        dataPathOut   % directory for the data (do not end with \)
        jobNameOut    % will read structure from ampac out file  jobname.out

        % and store indo results in jobname.ido
        % Output success/fail
        indo_succeed;   % TRUE means Indo calcs succeeded, FALSE means an exception was thrown
        indo_err_msg;   % Error msg from output
        % Atomic basis set:
        norb       % number of atomic basis functions, and hf orbitals
        aorbAtom   % (1,i) atom on which the ith atomic orbital resides
        aorbType   % (1,i) type of ith atomic orbital
        %  {s1=0,s2=1,p2x=2,p2y=3,p2z=4,s3=5,p3x=6,p3y=7,p3z=8}
        % Hartree Fock Results:
        nfilled    % number of filled molecular orbitals
        hfE        % HF ground state energy
        orbE       % (1,i) energy of ith orbital
        orb        % (i,j) ith component of jth orbital
        indoOutput % output from indo
        % SCI results
        nsci        % number of sci states, with first being ground state
        nscibasis   % number of basis functions (first being ground state)
        esci        % (1,i) energies of ith sci state
        r           %( i,j, icomp) transition (position) operator icomp = x,y,z
        wfsci       % (i,j)  ith component of jth state
        ehsci       % (i,1) hole of the ith SCI basis function (0 if GS)
        % (i,2) = elec of the ith SCI basis function (0 if GS)

        indoExe;    % Path to INDO exe
    end % properties
    properties (Transient)
        osc         % (1,i) oscillator strength from gs to state i
        rx          % r(:,:,1) for backwards compatibility, don't use now
        ry          % r(:,:,2)
        rz          % r(:,:,3)
    end 
    methods (Access = private)
        function readIndo(obj, inputfile)
            if (nargin > 1)
                filename = inputfile;
            else
                filename = [obj.dataPathOut,'\',obj.jobNameOut,'.ido'];
            end
            
            fid1 = fopen(filename);
            if (fid1 == -1)
               error(['in Indo.readIndo, could not find file: ',filename]);
            end
            obj.norb = fread(fid1,1,'integer*4');
            obj.aorbAtom = fread(fid1,[1,obj.norb],'integer*4');
            obj.aorbAtom = obj.aorbAtom +1; % C++ starts count at 0 instead of 1
            obj.aorbType = fread(fid1,[1,obj.norb],'integer*4');

            ntest = fread(fid1,1,'integer*4');
            if (ntest ~= obj.norb)
               error('atomic and fock basis sizes differ');
            end
            obj.nfilled = fread(fid1,1,'integer*4');
            obj.hfE  = fread(fid1,1,'real*8');
            obj.orbE = fread(fid1,[1,obj.norb],'real*8');
            obj.orb = fread(fid1,[obj.norb,obj.norb],'real*8');

            obj.nsci = fread(fid1,1,'integer*4');
            obj.nscibasis = fread(fid1,1,'integer*4');
            obj.esci = fread(fid1,[1,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r = zeros(obj.nsci,obj.nsci,3);
            obj.r(:,:,1) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r(:,:,2) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r(:,:,3) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            temp = fread(fid1,[2,obj.nscibasis],'integer*4');
            obj.ehsci = temp' +1; % +1 fixes the counting from 0 issue
            obj.wfsci = fread(fid1,[obj.nscibasis,obj.nsci],'real*8');
            fclose(fid1);
        end
    end
    methods (Static)
        function res = defaultConfig()
            res.charge = 1;
            res.norbs = 100;
            res.nstates = 25;
            res.field = [0,0,0];
            res.initial_shiftc = 0.0;
            res.initial_shift_step = 0.0;
            res.min_shift_step = 0.0;
            res.max_shift_step = 0.0;
            res.initial_second_shift_step = 0.0;
            res.min_second_shift_step = 0.0;
            res.max_second_shift_step = 0.0;
            res.initial_eeint = 1.0;
            res.initial_eestep = 0.0;
            res.min_eestep = 0.0;
            res.max_eestep = 0.0;
            res.initial_conv = 1e-10;
            res.min_conv = 1e-10;
            res.max_inner_iter = 300;
            res.max_iter = 10000;
            res.dm_guess = 'default';
            res.try_default_first = false;
            res.output_dm = false;
            res.pot_file = [];
        end
        
        function obj = LoadExistingData(filename, varargin)
            try
            if (~exist(filename,'file'))
                throw(MException('Indo:LoadExistingData:FileDNE',...
                    'Indo: Cannot locate INDO file'));
            end
            catch exception
                disp(exception);
            end
            
            obj = Indo();
            
            listing = dir(filename);
            if (listing.bytes == 0)
                obj.indo_succeed = false;
                obj.indo_err_msg = 'File size is 0 bytes. Calc failed on previous run.';
            else
                obj.readIndo(filename);
                obj.indo_succeed = true;
                obj.indo_err_msg = [];
            end         
            
            if (nargin < 2)
                obj.config = [];
            else
                obj.config = varargin{1};
            end
            if (nargin < 3)
                obj.jobNameIn = [];
            else
                obj.jobNameIn = varargin{2};
            end
            if (nargin < 4)
                obj.dataPathIn = [];
            else
                obj.dataPathIn = varargin{3};
            end
            if (nargin < 5)
                obj.jobNameOut = [];
            else
                obj.jobNameOut = varargin{4};
            end
            if (nargin < 6)
                obj.dataPathOut = [];
            else
                obj.dataPathOut = varargin{5};
            end
        end
    end
    methods       
        function res = Indo(ConfigIn, dataPathIn, jobNameIn, dataPathOut, jobNameOut, sysvars)
            
            if (nargin < 1)
                res.config = Indo.defaultConfig();
            else
                res.config = ConfigIn;
            end
            if (nargin < 2)
                res.dataPathIn = 'data';
            else
                res.dataPathIn = dataPathIn;
            end
            if (nargin < 3)
                res.jobNameIn = 'jobname';
            else
                res.jobNameIn = jobNameIn;
            end
            if (nargin < 4)
                res.dataPathOut = res.dataPathIn;
            else
                res.dataPathOut = dataPathOut;
            end
            if (nargin < 5)
                res.jobNameOut = res.jobNameIn;
            else
                res.jobNameOut = jobNameOut;
            end
            if (nargin < 6)
                sysvars = ECESysVars.getInstance;
                warning('off','ECESysVars:AlreadyInit');
                sysvars.initialize;
                warning('on','ECESysVars:AlreadyInit');
            end
            res.indoExe = ['"',sysvars.getVars('indo'),'"'];
            
            if (nargin > 0)
%                 jobstring = [res.indoExe,' "',res.dataPath,'\',res.jobName,'"', ...
%                     ' ',num2str(res.config.charge), ...
%                     ' ',num2str(res.config.norbs), ...
%                     ' ',num2str(res.config.nstates)];
%                 if (norm(res.config.field) > 0.0)
%                     jobstring = [jobstring, ...
%                         ' ', num2str(res.config.field(1)), ...
%                         ' ', num2str(res.config.field(2)), ...
%                         ' ', num2str(res.config.field(3))];
%                 end
                
                % The following is a workaround until I can recompile INDO
                % to not barf with the new AMPAC file format which lists
                % the keyword "CARTESIAN" 3 times instead of 2 and INDO has
                % a tizzy
                workaround = 0;
                % This is to check to see if you're using an AMPAC file or Gaussian
                if (exist(fullfile(res.dataPathIn,[res.jobNameIn,'.out']),'file'))
                    workaround = 1;
                    out_file = fileread(fullfile(res.dataPathIn,[res.jobNameIn,'.out']));
                    out_file = textscan(out_file,'%s','delimiter','\n');
                    out_file = out_file{1};
                    for i = 1:length(out_file)
                        if (~isempty(regexpi(out_file{i},'AMPAC Version','match')))
                            vers = regexp(out_file{i},'Version (?<v>\d*)','names');
                            break;
                        end
                    end
                    vers = str2double(vers.v);
                    if (vers > 9.9)   % If AMPAC version 10 or higher
                        [~,tmp_name,~]= fileparts(tempname);
                        tmp_name = ['indo_', tmp_name];
                        while (exist(fullfile(res.dataPathIn, [tmp_name,'.out']),'file'))
                            [~,tmp_name,~] = fileparts(tempname);
                            tmp_name = ['indo_', tmp_name];
                        end
                        fid = fopen(fullfile(res.dataPathIn,[tmp_name,'.out']),'w');
                        for i = 1:length(out_file)
                            if (~isempty(regexp(out_file{i},...
                                    'CARTESIAN COORDINATES READ IN BUT CALCULATION TO BE RUN IN INTERNAL COORDINATES','match')))
                                fprintf(fid,'%s\n','');
                            else
                                fprintf(fid,'%s\n',out_file{i});
                            end
                        end
                        fclose(fid);
                        old_name = res.jobNameIn;
                        res.jobNameIn = tmp_name;
                    end
                end
                % END PART 1 WORKAROUND
                
                ipf_file = res.write_param_file();
                jobstring = [res.indoExe,' "',res.dataPathOut,'\',ipf_file,'"'];
                                
                % disp(['about to do: ',jobstring]);
                [~, result] = system(jobstring);
                % system(jobstring);  % For diagnostics, this will output everything
                
                % Check to see if we wrote anything to the IDO file...if
                % not we died somewhere
                pause(2);
                listing = dir(fullfile(res.dataPathOut,[res.jobNameOut,'.ido']));
                if (listing.bytes == 0)
                    disp(result);
                    throw(MException('INDO:IndoFailure',...
                        'INDO: The program has failed to output data. Above is the command window output'));
                end
                
                delete([res.dataPathOut,'\',ipf_file]);
                
                % PART 2 WORKAROUND
                if (workaround)
                    if (vers > 9.9)
                        res.jobNameIn = old_name;
                        delete(fullfile(res.dataPathIn, [tmp_name,'.out']));
                    end
                end
                % END PART 2 WORKAROUND
                
                res.indoOutput = result;
                if (~isempty(regexpi(result, '.*Repulse exception.*', 'start')))
                    [tok, ~] = regexpi(result, '.*Repulse exception(.*)', 'tokens');
                    res.indo_err_msg = num2str(cell2mat(tok{1}(1)));
                    res.indo_succeed = false;
                    res.norb = [];
                    res.aorbAtom = [];
                    res.aorbType = [];
                    res.nfilled = [];
                    res.hfE = [];
                    res.orbE = [];
                    res.orb = [];
                    res.nsci = [];
                    res.nscibasis = [];
                    res.esci = [];
                    res.r = [];
                    res.wfsci = [];
                    res.ehsci = [];
                else
                    res.indo_succeed = true;
                    res.readIndo();
                end
            end
        end % INDO constructor
        function res = get_osc(obj)
            res = zeros(1,obj.nsci);
            for i=1:obj.nsci
                res(1,i) = (obj.esci(1,i)-obj.esci(1,1)) * ...
                   ( obj.r(1,i,1)^2 + obj.r(1,i,2)^2 + obj.r(1,i,3)^2 );
            end
        end
        function res = get_rx(obj)
            res = obj.r(:,:,1);
        end
        function res = get_ry(obj)
            res = obj.r(:,:,2);
        end
        function res = get_rz(obj)
            res = obj.r(:,:,3);
        end
        function res = dipole(obj,istate,jstate)
            % returns a vector that is the dipole(istate=jstate)
            % or transition moment (istate ~= jstate)
            res = reshape(obj.r(istate,jstate,:),[3,1]);
        end
        function ipf_name = write_param_file(obj)
            [~,ipf_name,~]= fileparts(tempname);
            ipf_name = ['indo_', ipf_name, '.ipf'];
            while (exist([obj.dataPathOut, '\', ipf_name],'file'))
                [~,ipf_name,~]= fileparts(tempname);
                ipf_name = ['indo_', ipf_name, '.ipf'];
            end
            fid = fopen([obj.dataPathOut, ipf_name], 'w');
            
            if (fid == -1)
                throw(MException('Indo:Write_Param_File:IOError','Indo: Cannot create param file.'));
            end
            
            fprintf(fid, '%s\r\n', ['jobname_in = "', obj.dataPathIn, '\', obj.jobNameIn, '"']);
            fprintf(fid, '%s\r\n', ['jobname_out = "', obj.dataPathOut, '\', obj.jobNameOut, '"']);

            if (strcmpi(obj.config.dm_guess, 'default'))
                fprintf(fid, '%s\r\n', 'dm_guess = default');
            else
                fprintf(fid, '%s\r\n', ['dm_guess = "', obj.config.dm_guess, '"']);
            end
            
            if (~isempty(obj.config.pot_file))
                fprintf(fid, '%s\r\n', ['pot_file = "', obj.config.pot_file, '"']);
            end
            
            fprintf(fid, '%s%i\r\n', 'charge = ', obj.config.charge);
            fprintf(fid, '%s%u\r\n', 'norbs = ', obj.config.norbs);
            fprintf(fid, '%s%u\r\n', 'nstates = ', obj.config.nstates);
            fprintf(fid, '%s%u\r\n', 'max_inner_iter = ', obj.config.max_inner_iter);
            fprintf(fid, '%s%u\r\n', 'max_iter = ', obj.config.max_iter);
            fprintf(fid, '%s%u\r\n', 'try_default_first = ', obj.config.try_default_first);
            fprintf(fid, '%s%u\r\n', 'output_dm = ', obj.config.output_dm);
            
            fprintf(fid, '%s%6.10f\r\n', 'efieldx = ', obj.config.field(1));
            fprintf(fid, '%s%6.10f\r\n', 'efieldy = ', obj.config.field(2));
            fprintf(fid, '%s%6.10f\r\n', 'efieldz = ', obj.config.field(3));
            
            fprintf(fid, '%s%3.10f\r\n', 'initial_shiftC = ', obj.config.initial_shiftc);
            fprintf(fid, '%s%3.10f\r\n', 'initial_shift_step = ', obj.config.initial_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'min_shift_step = ', obj.config.min_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'max_shift_step = ', obj.config.max_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'initial_second_shift_step = ', obj.config.initial_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'min_second_shift_step = ', obj.config.min_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'max_second_shift_step = ', obj.config.max_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'initial_eeint = ', obj.config.initial_eeint);
            fprintf(fid, '%s%3.10f\r\n', 'initial_eestep = ', obj.config.initial_eestep);
            fprintf(fid, '%s%3.10f\r\n', 'max_eestep = ', obj.config.max_eestep);
            fprintf(fid, '%s%3.10f\r\n', 'min_eestep = ', obj.config.min_eestep);
            
            fprintf(fid, '%s%3.10e\r\n', 'initial_conv = ', obj.config.initial_conv);
            fprintf(fid, '%s%3.10e\r\n', 'min_conv = ', obj.config.min_conv);
            
            fclose(fid);
        end
            
    end % methods
end % class