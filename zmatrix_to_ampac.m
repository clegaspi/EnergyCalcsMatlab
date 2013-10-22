function [ succeed ] = zmatrix_to_ampac(zmatrix, filepath, filename_noext, varargin)
%ZMATRIX_TO_AMPAC Summary of this function goes here
%   Detailed explanation goes here
    natoms = size(zmatrix, 1);
    
    if (nargin > 3)
        overwrite = varargin{1};
    else
        overwrite = true;
    end
    
    if (nargin > 4)
        params = varargin{2};
    else
        params = 'AM1 rhf singlet 1scf t=auto geo=ok';  % Single point calculation
    end
    
    if (exist(fullfile(filepath,[filename_noext,'.dat']), 'file') && ~overwrite)
        succeed = 2;
        return;
    end
    
    try
        fid = fopen(fullfile(filepath,[filename_noext,'.dat']),'w');

        fprintf(fid, '%s\r\n%s\r\n%s\r\n', params, 'Title', 'Comment');
        % fprintf(fid, '%s\r\n\r\n\r\n', 'INSERT KEYWORDS HERE');
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

        fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s', '0', '0.000000', '0', '0.000000', '0', '0.000000', '0', '0','0','0');

        fclose(fid);
    catch exception
        succeed = 0;
        return;
    end
    
    succeed = 1;
end