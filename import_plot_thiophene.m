path = 'C:\Users\clegaspi\Documents\Thiophenes\INDO\';
mols = {'2T-opt', '2T-planar-opt', ...
    '2T-33-dichloro-opt', '2T-33-dichloro-planar-opt', ...
    '2T-33-dimethyl-opt', '2T-33-dimethyl-planar-opt', ...
    '2T-44-dichloro-opt', '2T-44-dichloro-planar-opt', ...
    '2T-44-dimethyl-opt', '2T-44-dimethyl-planar-opt'};

mis_exen = zeros(1,length(mols));
gs_en = zeros(1,length(mols));

for i = 1:length(mols)
    indo = Indo.LoadExistingData([path,mols{i},'.ido']);
    % eexc = indo.esci - indo.esci(1);
    eexc = indo.esci;
    gs_en(i) = indo.esci(1);
    osc = indo.get_osc();
    figure(1);
    hold on;
    plot(i*ones(1,length(eexc)), eexc, 'g^');
    for j = 1:length(eexc)
        plot(i, eexc(j), 'bo', 'MarkerSize', 30*osc(j)/max(osc) + 1e-3);
    end
    idx = find((30*osc./max(osc) + 1e-3) > 10);
    if (isempty(idx))
        keyboard;
    end
    
    mis_exen(i) = eexc(idx(1));
end