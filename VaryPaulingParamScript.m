k=1;
S = load('s:\NoAngle\8-merPPV-0VA.mat');
myexp = S.obj;
myexp.data(1).load_to_memory('indo','load');
myexp.data(1).generate_ampac_file('out','p:\');

ppdata = repmat(struct(),4,1);
boandr = cell(4,1);
pp = [0.1, 0.2, 0.4, 0.5];
for i = 1:4
    [rgs, rex, deltar, ~, numruns, eszmat, ampac_energy] = OptExcStateStructure('p:\', myexp.data(1).ampac_hash,...
        'o:\', myexp.data(1).indo_hash,...
        'field', [0 0 0],...
        'indo', myexp.data(1).raw_indo,...
        'c', pp(i));
    boandr{i} = {rgs, rex, deltar, numruns, eszmat, ampac_energy};
    
    ppdata(i).exdeltaAM1{1,1}(k) = ampac_energy - myexp.data(1).Ehf;

    myexp.data(1).load_to_memory('ampac','load');
    amp = myexp.data(1).raw_ampac;
    amp.ampac_succeed = true;

    ppdata(i).exindo = Indo.LoadExistingData(['o:\', myexp.data(1).indo_hash, '-new.ido'],[],[],[]);

    tmp = ppdata(i).exindo.dipole(1,1);
    ppdata(i).exdpgs{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = ppdata(i).exindo.dipole(2,2);
    ppdata(i).exdp{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    tmp = ppdata(i).exindo.dipole(3,3);
    ppdata(i).exdp2exc{1,1}(k) = sum(tmp .^ 2, 1) .^ (0.5);

    exnfill = ppdata(i).exindo.nfilled;
    tmp = ppdata(i).exindo.orbE([exnfill exnfill+1]);
    ppdata(i).exhlgap{1,1}(k) = tmp(2) - tmp(1);

    tempds = ECEDataStruct(amp, ppdata(i).exindo,[], [], [], [], []);
    ppdata(i).exeexc{1,1}(k,1:25) = tempds.Eexc(:);

    ppdata(i).exigs{1,1}(k) = ppdata(i).exindo.esci(1);

    ppdata(i).exopint{1,1}(k,1:25) = tempds.Tint(:);
    tmp = squeeze(ppdata(i).exindo.r(2,:,:));
    for l = 1:25
        ppdata(i).extpint{1,1}(k,l) = (ppdata(i).exeexc{1,1}(k,l) - ppdata(i).exeexc{1,1}(k,2)) * sum(tmp(l,:) .^ 2);
    end
end

%% Plot
figure(8)
hold on

ppidx = 4;

for ppidx = 1:4
    xaxis = pp(ppidx);
    
    % GS Stuff
    
    for k = 1:length(xaxis)
        plot(xaxis(k), eexc{1,1}(k,:) + igs{1,1}(k), 'g^');
    end

    maxopint = max(opint{1,1}(:));
    for k = 1:length(xaxis)
        for l = 1:25
            plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k), 'bo', 'MarkerSize', (opint{1,1}(k,l)*30 / maxopint) + 1e-3);
        end
    end

    maxtpint = max(tpint{1,1}(:));

    for k = 1:length(xaxis)
        for l = 1:25
            if (tpint{1,1}(k,l) > 0)
                plot(xaxis(k), eexc{1,1}(k,l) + igs{1,1}(k), 'rs', 'MarkerSize', (tpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
            end
        end
    end

    % Excited state stuff

    for k = 1:length(xaxis)
        plot(xaxis(k), ppdata(ppidx).exeexc{1,1}(k,:) + igs{1,1}(k) + ppdata(ppidx).exdeltaAM1{1,1}(k), 'k^');
    end

    % maxexopint = max(exopint{1,1}(:));
    for k = 1:length(xaxis)
        for l = 1:25
            plot(xaxis(k), ppdata(ppidx).exeexc{1,1}(k,l) + igs{1,1}(k) + ppdata(ppidx).exdeltaAM1{1,1}(k), 'mo', 'MarkerSize', (exopint{1,1}(k,l)*30 / maxopint) + 1e-3);
        end
    end

    % maxtpint = max(tpint{1,1}(:));

    for k = 1:length(xaxis)
        for l = 1:25
            if (tpint{1,1}(k,l) > 0)
                plot(xaxis(k), ppdata(ppidx).exeexc{1,1}(k,l) + igs{1,1}(k) + ppdata(ppidx).exdeltaAM1{1,1}(k), 'cs', 'MarkerSize', (extpint{1,1}(k,l)*30 / maxtpint) + 1e-3);
            end
        end
    end
end