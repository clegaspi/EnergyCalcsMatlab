function successful = planarize_zmatrix(in_file, out_file)
    % This takes the zmatrix of a molecule that you want the carbons to be
    % planar and planarizes them. You should open in AMPAC and delete the
    % hydrogens and add them again before processing.
    successful = 0;
    zmat = ampac_to_zmatrix(in_file);
    for i = 1:size(zmat,1)
        if (strcmpi(zmat{i,2},'C'))
            if (abs(str2double(zmat{i,5})) > 160 && abs(str2double(zmat{i,5})) < 181)
                zmat{i,5} = '180.000000';
            else
                zmat{i,5} = '0.000000';
            end
        end
    end
    [out_path, out_fn, ~] = fileparts(out_file);
    successful = zmatrix_to_ampac(zmat,out_path,out_fn);
    if (successful)
        disp('Be sure to open the file and delete the hydrogens and add them again for good placement.');
    end
end