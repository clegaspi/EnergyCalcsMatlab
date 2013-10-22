classdef StateCoupling < handle
    %STATECOUPLING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lower_state;
        upper_state;
        upper_bound;
        lower_bound;
        coupling;
        field_coupling;
        efv;
        jobpath_in;
        jobname_in;
        jobpath_out;
        jobname_out;
        final_indo;
    end
    
    methods
        function obj = StateCoupling(lower_state,upper_state,lb,ub,field_vector,jobpath_in,jobname_in,varargin)
            if (lb == ub)
                throw(MException('StateCoupling:InvalidBounds','StateCoupling: Bounds must not be the same!'));
            end
            if (lb > ub)
                obj.upper_bound = lb;
                obj.lower_bound = ub;
            else
                obj.upper_bound = ub;
                obj.lower_bound = lb;
            end
            if (lower_state == upper_state)
                throw(MException('StateCoupling:InvalidStates','StateCoupling: States must not be the same!'));
            end
            if (lower_state > upper_state)
                obj.upper_state = lower_state;
                obj.lower_state = upper_state;
            else
                obj.upper_state = upper_state;
                obj.lower_state = lower_state;
            end
            obj.efv = field_vector ./ norm(field_vector);
            obj.jobpath_in = jobpath_in;
            obj.jobname_in = jobname_in;
            if (length(varargin) == 2)
                % This is to specify the output file name for the INDO
                % file for the coupling calculation. Necessary for parallel
                % calculation!
                obj.jobpath_out = varargin{1};
                obj.jobname_out = ['coup_',varargin{2}];
            else
                obj.jobpath_out = jobpath_in;
                obj.jobname_out = ['coup_',jobname_in];
            end
        end
        
        function [dE, fs] = run(obj)
            options = optimset('Display','iter','TolX',1e-8);
            [fs, dE, exitflag] = fminbnd(@obj.get_coupling, obj.lower_bound, obj.upper_bound, options);
            obj.coupling = dE;
            obj.field_coupling = fs;
        end
        
        function dE = get_coupling(obj, fs)
            res = [];
            res.charge = 0;
            res.norbs = 500;
            res.nstates = 25;
            res.field = fs .* obj.efv;
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
            
            if (exist([obj.jobpath_in, obj.jobname_in, '-dm.bin'],'file'))
                res.dm_guess = [obj.jobpath_in, obj.jobname_in, '-dm.bin'];
            else
                res.dm_guess = 'default';
            end
            
            res.try_default_first = true;
            res.output_dm = true;
            res.pot_file = [];

            obj.final_indo = Indo(res, obj.jobpath_in, obj.jobname_in, obj.jobpath_out, obj.jobname_out);
            
            dE = abs(obj.final_indo.esci(obj.upper_state)-obj.final_indo.esci(obj.lower_state));
        end   
        
            
    end
    
end

