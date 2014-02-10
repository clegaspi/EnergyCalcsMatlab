function [ varargout ] = ECE_system_vars( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    global INDOEXE;
    global AMPACEXE;
    global DEFAULT_POOL_SIZE;
    global ECE_DATA_PATH;
    output_args = [];
    
    getvars = ~isempty(find(cellfun(@(x)strcmpi(x,'get'), varargin),1)) || nargin == 0;
    checkvars = ~isempty(find(cellfun(@(x)strcmpi(x,'check'), varargin),1));
    
    if (isempty(getenv('F_EM64T_REDIST11')) && isempty(getenv('F_IA32_REDIST11')) && isempty(getenv('F_IA64_REDIST11')))
        if (checkvars)
            output_args.fortran = 0;
        else
            throw(MException('ECE_SYSTEM_VARS:NoFortran','ECE_SYSTEM_VARS: Fortran Libraries not found!'));
        end
    else
        output_args.fortran = 1;
    end
        
    if (~getvars)
        setall = ~isempty(find(cellfun(@(x)strcmpi(x,'setall'), varargin),1));
        setrest = ~isempty(find(cellfun(@(x)strcmpi(x,'setrest'), varargin),1));
        
        indoidx = find(cellfun(@(x)strcmpi(x,'indo'), varargin));
        ampacidx = find(cellfun(@(x)strcmpi(x,'ampac'), varargin));
        poolidx = find(cellfun(@(x)strcmpi(x,'poolsize'), varargin));
        dataidx = find(cellfun(@(x)strcmpi(x,'datapath'), varargin));
        
        set_indo = setall || (~isempty(indoidx) && strcmpi(varargin{indoidx+1},'set'))...
            || (isempty(indoidx) && setrest);
        set_ampac = setall || (~isempty(ampacidx) && strcmpi(varargin{ampacidx+1},'set'))...
            || (isempty(ampacidx) && setrest);
        set_pool = setall || (~isempty(poolidx) && strcmpi(varargin{poolidx+1},'set'))...
            || (isempty(poolidx) && setrest);
        set_data = setall || (~isempty(dataidx) && strcmpi(varargin{dataidx+1},'set'))...
            || (isempty(dataidx) && setrest);

        if (set_indo)
            if (exist('c:\mscpp\demo-dci\Release\demo-dci.exe','file'))
                INDOEXE = 'c:\mscpp\demo-dci\Release\demo-dci.exe';
            else
                throw(MException('ECE_SYSTEM_VARS:NoINDO','ECE_SYSTEM_VARS: INDO not found!'));
            end
        elseif (checkvars)
            if (~exist(INDOEXE,'file'))
                output_args.indo = [];
            end
        elseif (~isempty(indoidx))
            if (exist(varargin{indoidx+1},'file'))
                INDOEXE = varargin{indoidx+1};
            else
                throw(MException('ECE_SYSTEM_VARS:NoINDO','ECE_SYSTEM_VARS: INDO not found!'));
            end
        end
        
        if (set_ampac)
            if (exist('c:\Program Files\Semichem, Inc.\Ampac-10.1\ampac.exe','file'))
                AMPACEXE = 'c:\Program Files\Semichem, Inc.\Ampac-10.1\ampac.exe';
            elseif (exist('c:\Program Files (x86)\Semichem, Inc.\Ampac-10.1\ampac.exe','file'))
                AMPACEXE = 'c:\Program Files (x86)\Semichem, Inc.\Ampac-10.1\ampac.exe';
            else
                throw(MException('ECE_SYSTEM_VARS:NoAmpac','ECE_SYSTEM_VARS: AMPAC not found!'));
            end
        elseif (checkvars)
            if (~exist(AMPACEXE,'file'))
                output_args.ampac = [];
            end
        elseif (~isempty(ampacidx))
            if (exist(varargin{ampacidx+1},'file'))
                INDOEXE = varargin{ampacidx+1};
            else
                throw(MException('ECE_SYSTEM_VARS:NoAmpac','ECE_SYSTEM_VARS: AMPAC not found!'));
            end
        end
        
        if (set_pool)
            DEFAULT_POOL_SIZE = str2double(getenv('NUMBER_OF_PROCESSORS')) - 1;
        elseif (checkvars)
            if (DEFAULT_POOL_SIZE > str2double(getenv('NUMBER_OF_PROCESSORS')))
                output_args.poolsize = [];
            end
        elseif (~isempty(poolidx))
            nproc = str2double(getenv('NUMBER_OF_PROCESSORS'));
            if (isnumeric(varargin{poolidx+1}) && varargin{poolidx+1} <= nproc && varargin{poolidx+1} > 0)
                DEFAULT_POOL_SIZE = int32(varargin{poolidx+1});
            else
                throw(MException('ECE_SYSTEM_VARS:InvalidPool','ECE_SYSTEM_VARS: Invalid pool size!'));
            end
        end
        
        if (set_data)
            if (exist(fullfile(pwd, '..\data'),'dir'))
                ECE_DATA_PATH = fullfile(pwd, '..\data');
            elseif (exist(fullfile(pwd, 'data'),'dir'))
                ECE_DATA_PATH = fullfile(pwd, 'data');
            else
                throw(MException('ECE_SYSTEM_VARS:InvalidDataPath','ECE_SYSTEM_VARS: Invalid data path!'));
            end
        elseif (checkvars)
            if (~exist(ECE_DATA_PATH,'dir'))
                output_args.datapath = [];
            end
        elseif (~isempty(dataidx))
            if (exist(varargin{dataidx+1},'dir'))
                DEFAULT_POOL_SIZE = fullpath(varargin{dataidx+1});
            else
                throw(MException('ECE_SYSTEM_VARS:InvalidDataPath','ECE_SYSTEM_VARS: Invalid data path!'));
            end
        end
    end
    
    if (~isfield(output_args,'indo'))
        output_args.indo = INDOEXE;
    end
    if (~isfield(output_args,'ampac'))
        output_args.ampac = AMPACEXE;
    end
    if (~isfield(output_args,'poolsize'))
        output_args.poolsize = DEFAULT_POOL_SIZE;
    end
    if (~isfield(output_args,'datapath'))
        output_args.datapath = ECE_DATA_PATH;
    end
    
    if (checkvars)
        if (~output_args.fortran || isempty(output_args.indo) || isempty(output_args.ampac)...
                || isempty(output_args.poolsize) || isempty(output_args.datapath))
            varargout{1} = false;
        else
            varargout{1} = true;
        end
        if (nargout == 2)
            varargout{2} = output_args;
        end
    else
        varargout{1} = output_args;
    end
end

