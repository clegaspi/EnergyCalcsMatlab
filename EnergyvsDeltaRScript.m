S = load('s:\NoAngle\8-merPPV-0.101VA.mat');
myexp = S.obj;
thisdat = myexp.data(1);
thisdat.load_to_memory('indo','load');
thisdat.generate_ampac_file('out','p:\');

[~,~,~, eszmat, ~, ~, ring_info, ~] = OptExcStateStructure('p:\', thisdat.ampac_hash,...
    'o:\', thisdat.indo_hash,...
    'field', [0 0 0],...
    'indo', thisdat.raw_indo,...
    'readifexist',...
    'algorithm','paulingwithAM1');

atoms_in_bonds = get_connectivity(eszmat);
[gs_xyz, ~] = ampac_to_xyz(['p:\',thisdat.ampac_hash,'.out']);
[ex_xyz, ~] = ampac_to_xyz(eszmat);

oldblsbyatom = zeros(size(eszmat,1));
deltarbyatom = zeros(size(eszmat,1));

for i = 1:length(atoms_in_bonds)
    for j = atoms_in_bonds{i}
        oldbl = sum((gs_xyz(i,:) - gs_xyz(j,:)) .^ 2) ^ 0.5;
        oldblsbyatom(i,j) = oldbl;
        oldblsbyatom(j,i) = oldbl;

        dr = (sum((ex_xyz(i,:) - ex_xyz(j,:)) .^ 2) ^ 0.5) - oldbl;
        deltarbyatom(i,j) = dr;
        deltarbyatom(j,i) = dr;
    end
end

res = [];
res.charge = 0;
res.norbs = 500;
res.nstates = 25;
res.field = [0 0 0];
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
res.dm_guess = ['o:\',thisdat.indo_hash, '-dm.bin'];
res.try_default_first = true;
res.output_dm = true;
res.pot_file = [];
            
deltar_const = 0:0.1:2;
Egs = zeros(1,length(deltar_const));
dEindo = zeros(1,length(deltar_const));

for i = 1:length(deltar_const)
    newblsbyatom = oldblsbyatom + deltar_const(i) * deltarbyatom;
    newzmatrix = distort_geometry(eszmat,ring_info,newblsbyatom,'am1optimize');
    zmatrix_to_ampac(newzmatrix,'p:\',[thisdat.indo_hash,'-deltarvar']);
    [~,~] = system(['"c:\Program Files\Semichem, Inc\Ampac-9.2\ampac.exe" "p:\',thisdat.indo_hash,'-deltarvar.dat"']);
    ampac = parseAmpac(['p:\',thisdat.indo_hash,'-deltarvar']);
    Egs(i) = ampac.Hf / 23.05;
    newindo = Indo(res, 'p:\', [thisdat.indo_hash,'-deltarvar']);
    dEindo(i) = newindo.esci(1,2) - newindo.esci(1,1);
end

%% Plot

figure(6);
hold on;
plot(deltar_const, Egs, '--k');
plot(deltar_const, dEindo, '--b');
plot(deltar_const, Egs + dEindo, '--r');


%% Calculate polaron relaxation energy

fes = Egs + dEindo;
[minval,minidx] = min(fes);
Epol = fes(1) - minval;
disp(['The polaron relaxation energy is ', num2str(Epol), ' eV with delta(R) const = ', ...
    num2str(deltar_const(minidx))]);
    
    

