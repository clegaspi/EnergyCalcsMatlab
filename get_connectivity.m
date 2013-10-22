function [ varargout ] = get_connectivity( input )
%GET_CONNECTIVITY Summary of this function goes here
%   Detailed explanation goes here
    
    OpenBabelEXE = 'C:\Program Files (x86)\OpenBabel-2.3.2\babel.exe';
    
    if (~exist(OpenBabelEXE, 'file'))
        throw(MException('get_connectivity:OpenBabelNotFound', 'get_connectivity: OpenBabel EXE not found. Please edit.'));
    end
    
    if (iscell(input))
        zmatrix = input;
        full_path = [cd, '\temp', num2str(int32(rand*1000))];
    else
        zmatrix = ampac_to_zmatrix(input);
        full_path = input(1:end-4);
    end
    
    natoms = size(zmatrix, 1);
    
    fid = fopen([full_path, '-get_conn.mopin'],'w');

    fprintf(fid, '%s\r\n\r\n\r\n', 'INSERT KEYWORDS HERE');
    % fprintf(fid, '%s\r\n\r\n', 'Title');
    % '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n'
    fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{1,2}, '0.000000', '1', '0.000000', '1', '0.000000', '1', '0','0','0');

    if (natoms > 1)
        fprintf(fid,'%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{2,2}, zmatrix{2,3}, '1',...
            '0.000000', '1', '0.000000', '1', zmatrix{2,6},'0','0');
    end
    if (natoms > 2)
        fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{3,2}, zmatrix{3,3}, '1',...
            zmatrix{3,4}, '1', '0.000000', '1', zmatrix{3,6}, zmatrix{3,7},'0');
    end
    if (natoms > 3)
        for i = 4:natoms
            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{i,2}, zmatrix{i,3}, '1',...
                zmatrix{i,4}, '1', zmatrix{i,5}, '1', zmatrix{i,6}, zmatrix{i,7}, zmatrix{i,8});
        end
    end

    fclose(fid);

    [~, ~] = system(['"', OpenBabelEXE, '" -imopin "', full_path, '-get_conn.mopin" -omol2 "',...
        full_path, '-get_conn.mol2"']);

    mf = fileread([full_path, '-get_conn.mol2']);

    mf = textscan(mf,'%s','delimiter','\n');
    mf = mf{1};

    bidx = find(cellfun(@(x)strcmp(x, '@<TRIPOS>BOND'), mf)) + 1;

    mf = {mf{bidx:end}};

    mf = cellfun(@(x)textscan(x,'%s'), mf);

    atomsbonds = cell(natoms,1);

    for i = 1:length(mf)
        atomsbonds{str2double(mf{i}(2))}(end+1) = str2double(mf{i}(3));
        atomsbonds{str2double(mf{i}(3))}(end+1) = str2double(mf{i}(2));
    end
    
    delete([full_path, '-get_conn.mol2']);
    delete([full_path, '-get_conn.mopin']);
    
    varargout{1} = atomsbonds;
    if (nargout > 1)
        varargout{2} = mf;
    end

end

