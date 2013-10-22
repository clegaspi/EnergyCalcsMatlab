function [ dipole ] = par_get_dipole( fs, myexp, nstates, all_data )
%GET_DIPOLE Summary of this function goes here
%   Detailed explanation goes here
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
    
    res.field = fs .* myexp.params.Efield_vector;
    randstr = strtrim(sprintf('%3.0f',rand*1000));
    
    copyfile([myexp.ampacdatapath, myexp.data(1).ampac_hash, '.out'],...
        [myexp.ampacdatapath, myexp.data(1).ampac_hash,'-',randstr,'.out']);
    myindo = Indo(res, myexp.ampacdatapath, [myexp.data(1).ampac_hash,'-',randstr],...
        myexp.ampacdatapath, [myexp.data(1).ampac_hash,'-',randstr]);

    delete([myexp.ampacdatapath, myexp.data(1).ampac_hash,'-',randstr,'.out']);
    delete([myexp.ampacdatapath, myexp.data(1).ampac_hash,'-',randstr,'.ido']);
    
    plot(repmat(fs,1,nstates), myindo.esci, 'c^')
    
    cross_num = find(arrayfun(@(x)dblcmp(x,fs),all_data));
        
    dipole = sum(myindo.r(cross_num+1,cross_num+1,:).^2).^0.5;

end

