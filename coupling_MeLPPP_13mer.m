%%
% Skip to plot? Set to false
if (true)
% Coupling parameters
path_to_nofield_matfile = [pwd,'\data\MeLPPP-13mer-N2Opt\Exp\NoAngle\MeLPPP-13mer-N2Opt-0VA.mat'];
path_coup_output = [pwd,'\data\MeLPPP-13mer-N2Opt\coupling\'];
fn_coupling = 'MeLPPP-13mer-N2Opt';
efield_vector_atoms = [1 186];
coup_par = {[2,3,0.038,0.039],...
    [3,4,0.043,0.044],...
    [4,5,0.044,0.045],...
    [5,6,0.048,0.049],...
    [6,7,0.0492,0.0493],...
    [7,8,0.0493,0.0494],...
    [8,9,0.0502,0.0504],...
    [9,10,0.0512,0.0516],...
    [10,11,0.0524,0.0526]};


if (matlabpool('size') == 0)
    matlabpool('3');
end

S = load(path_to_nofield_matfile);
myexp = S.obj;

myexp.data(1).load_to_memory('ampac','load');   % Load from ampac because the ampac_to_xyz gives diff cartesian?

efv = myexp.data(1).raw_ampac.r(:,efield_vector_atoms(1))-myexp.data(1).raw_ampac.r(:,efield_vector_atoms(2));

myexp.data(1).generate_ampac_file('out',path_coup_output);
movefile(fullfile(path_coup_output, [myexp.data(1).ampac_hash,'.out']),...
    fullfile(path_coup_output, [fn_coupling,'.out']));
if (exist(fullfile(myexp.data(1).indo_file_path, [myexp.data(1).indo_hash,'-dm.bin']), 'file'))
    copyfile(fullfile(myexp.data(1).indo_file_path, [myexp.data(1).indo_hash,'-dm.bin']),...
        fullfile(path_coup_output, [fn_coupling,'-dm.bin']));
end

coup = zeros(1,length(coup_par));
fs = coup;
indo = repmat(Indo(),1,length(coup_par));

%parfor idx = 1:length(coup_par)
parfor idx = 1:3
    disp(['Running Crossing #',num2str(2*idx)]);
    mycoup = StateCoupling(coup_par{2*idx}(1), coup_par{2*idx}(2),...
        coup_par{2*idx}(3), coup_par{2*idx}(4), efv,...
        path_coup_output,...
        fn_coupling,...
        path_coup_output,...
        [fn_coupling,'-Cross',num2str(2*idx)]);

    % [coup(idx), fs(idx)] = mycoup.run();
    [tc,tf] = mycoup.run()

    % indo(idx) = mycoup.final_indo;

    %disp(['Crossing ', num2str(idx)]);
    %disp(['Coupling: ', num2str(coup(idx))]);
    %disp(['Critical Field: ', num2str(fs(idx))]);
    
    disp(['Crossing ', num2str(2*idx)]);
    disp(['Coupling: ', num2str(tc)]);
    disp(['Critical Field: ', num2str(tf)]);
end

% save(fullfile(path_coup_output,[fn_coupling,'.mat']));

end

%% Plots
% Must have loaded state data already and plotted it in figure 8

ncross = length(coup_par);

figure(8);
hold on;

states_hdl = zeros(1,ncross);
op_hdl = cell(1,ncross);
ecenter = zeros(1,ncross);
mingap_energies = zeros(ncross,2);
mingap_hdl = zeros(1,ncross);

for i = 1:ncross
    nstates = indo(i).nsci;
    states_hdl(i) = plot(repmat(fs(i),1,nstates), indo(i).esci - igs{1,1}(1), 'm^');
    
    Tint = zeros(1,nstates);
    tmp = Tint;
    for j = 1:nstates
        Tint(j) = (indo(i).esci(j) - indo(i).esci(1)) * sum(indo(i).r(1,j,:).^2);
    end
    % maxopint should still be set
    for k = 1:nstates
        if (Tint(k)*30/maxopint > 1)
            tmp(k) = plot(fs(i), indo(i).esci(k) - igs{1,1}(1), 'bo', 'MarkerSize', (Tint(k)*30 / maxopint) + 1e-3);
        end
    end
    op_hdl{i} = tmp;

    ecenter(i) = ((indo(i).esci(coup_par{i}(1)) + indo(i).esci(coup_par{i}(2))) / 2) - igs{1,1}(1);
    mingap_energies(i,1) = indo(i).esci(coup_par{i}(1)) - igs{1,1}(1);
    mingap_energies(i,2) = indo(i).esci(coup_par{i}(2)) - igs{1,1}(1);
    mingap_hdl(i) = line(repmat(fs(i),1,2),mingap_energies(i,:),'Color','r','LineWidth',1.5);
end

coupling_hdl = plot(fs, ecenter, 'r*');

%% Delete plots

for i = 1:ncross
    delete(states_hdl(i));
    delete(mingap_hdl(i));
end

delete(coupling_hdl);