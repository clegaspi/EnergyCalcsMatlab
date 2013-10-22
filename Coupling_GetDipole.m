%% Get dipole moments after coupling points

% Coupling data AND energy data should be loaded before running this as
% well as running the cell which creates coupling and critfield vectors

exp_drive = 'j';
% exp_drive = 's';
mol_prefix = '17-merPPV-N2Optimized';
nstates = 25;

datapts = (critfield(1:end-1)+critfield(2:end))./2;
dipole = zeros(1,length(datapts)+1);
indo_save = repmat(Indo(),1,length(datapts));

% for i = 1:length(datapts)
%     datapts(i) = closest_member(datapts(i),field);
% end
% 
% [~,~,idx] = closest_member(critfield(end),field);
% 
% datapts = [datapts field(idx(1)+5)];
%     
% for i = 1:length(datapts)
%     S = load([exp_drive, ':\NoAngle\', mol_prefix, '-', num2str(datapts(i)), 'VA.mat']);
%     myexp = S.obj;
%      
%     myexp.data(1).load_to_memory('indo','load');
%     plot(repmat(datapts(i),1,nstates), myexp.data(1).raw_indo.esci, 'c^')
%     dipole(i) = sum(myexp.data(1).raw_indo.r(i+1,i+1,:).^2).^0.5;
% end

S = load([exp_drive, ':\NoAngle\', mol_prefix, '-',num2str(field(2)),'VA.mat']);
myexp = S.obj;

% numTasks = 6;
% [split, ~] = pctdemo_helper_split_vector(datapts, numTasks);
% 
% configName = defaultParallelConfig();
% sched = findResource('scheduler', 'Configuration', configName);
% job = createJob(sched);
% 
% for i = 1:numTasks
%     createTask(job, @par_get_dipole, 1, {split{i}, myexp, nstates, datapts});
% end
% 
% submit(job);
% waitForState(job, 'finished');
% y = getAllOutputArguments(job);
% destroy(job);
% cat(2, y{:})

    res = [];
    res.charge = 0;
    res.norbs = 1000;
    res.nstates = nstates;
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

    res.dm_guess = [myexp.indodatapath, myexp.data(1).indo_hash, '-dm.bin'];

    res.try_default_first = true;
    res.output_dm = false;
    res.pot_file = [];

for i = 1:length(datapts)
    res.field = datapts(i) .* myexp.params.Efield_vector;

%     randstr = strtrim(sprintf('%3.0f',rand*1000));
    
    myindo = Indo(res, myexp.ampacdatapath, myexp.data(1).ampac_hash);
    plot(repmat(datapts(i),1,nstates), myindo.esci, 'c^')
    drawnow;
    dipole(i) = sum(myindo.r(i+1,i+1,:).^2).^0.5;
    indo_save(i) = myindo;
end

% Need to add arbitrary point after last calculated crossing for that
% dipole

[~,~,idx] = closest_member(critfield(end),field);

datapts = [datapts field(idx(1)+1)];    % Adjust 5 to be after the coupling region and before any other coupling region

S = load([exp_drive, ':\NoAngle\', mol_prefix, '-', num2str(datapts(end)), 'VA.mat']);
myexp = S.obj;

myexp.data(1).load_to_memory('indo','load');
plot(repmat(datapts(end),1,nstates), myexp.data(1).raw_indo.esci, 'c^')
dipole(end) = sum(myexp.data(1).raw_indo.r(length(dipole)+1,length(dipole)+1,:).^2).^0.5;

save('C:\Users\Christian\Documents\Research\Yaron\dyes2\data\17-merPPV\Coupling\dipoles.mat',...
    'datapts', 'dipole', 'indo_save');


    
    
    
    