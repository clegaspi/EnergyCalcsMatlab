% mols = {'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-ortho.ido','r'; ...
%     'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-meta.ido','b';
%     'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-para.ido','k';
%     'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-UV-ortho.ido','y'; ...
%     'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-UV-meta.ido','c';
%     'C:\Users\clegaspi\Documents\Chlorine Test\diarylethene\DAE-UV-para.ido','g';};

mols = {'C:\Users\clegaspi\Documents\Thiophenes\INDO\2T-opt.ido','k'};

% clf;

for i = 1:size(mols,1)
    filepath = mols{i,1};
    color = mols{i,2};
    indo = Indo.LoadExistingData(filepath);

    exc = indo.esci-indo.esci(1);
    osc = indo.get_osc();
    prob = exc .* (osc.^2);

    % create a spectrum as a function of energy
    ener = 2.0:0.025:8.0; % eV
    spec = zeros(size(ener));
    LWother = 0.2; % assume a line width of 0.1eV for effects other than torsion
    figure(3)

    for it=1:length(prob)
       % we assume a gaussian centered at exc(i), with width LWother
       for ie = 1:length(ener)
          g = prob(it) * exp(-(ener-exc(it)).^2./LWother.^2);
          hold on;
          plot(1e7./(ener*8065.73) , g, color);
          spec = spec + g;
       end
    end
    % spec = spec./max(spec);
    figure(4)
    hold on
    plot(1e7./(ener*8065.73),spec,color);
end

figure(3);
xlabel('Wavelength (nm)');
ylabel('Absorbance (AU)');
title('Contributing Gaussians');
figure(4);
xlabel('Wavelength (nm)');
ylabel('Absorbance (AU)');
title('Summed Gaussians');