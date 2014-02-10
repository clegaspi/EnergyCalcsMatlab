function [ xyz, atom_type ] = ampac_to_xyz( input )
%AMPAC_TO_XYZ Summary of this function goes here
%   Detailed explanation goes here

    sysvars = ECESysVars.getInstance;
    ampacexe = sysvars.getVars('ampac');
    
    if (iscell(input))
        zmatrix = input;
    else
        zmatrix = ampac_to_zmatrix(input);
    end
    
    full_path = [pwd, '\temp', num2str(int32(rand*1000))];
    natoms = size(zmatrix,1);
    
    try
        params = 'am1 rhf singlet 1scf t=auto geo=ok';
        fid = fopen([full_path,'.dat'],'w');

        fprintf(fid, '%s\r\n%s\r\n%s\r\n', params, 'Title', 'Comment');

        fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{1,2}, '0.000000', '1', '0.000000', '1', '0.000000', '1', '0','0','0');

        if (natoms > 1)
            fprintf(fid,'%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{2,2}, zmatrix{2,3}, '0',...
                '0.000000', '1', '0.000000', '1', zmatrix{2,6},'0','0');
        end
        if (natoms > 2)
            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{3,2}, zmatrix{3,3}, '0',...
                zmatrix{3,4}, '1', '0.000000', '1', zmatrix{3,6}, zmatrix{3,7},'0');
        end
        if (natoms > 3)
            for i = 4:natoms
                fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{i,2}, zmatrix{i,3}, '0',...
                    zmatrix{i,4}, '1', zmatrix{i,5}, '1', zmatrix{i,6}, zmatrix{i,7}, zmatrix{i,8});
            end
        end

        fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s', '0', '0.000000', '0', '0.000000', '0', '0.000000', '0', '0','0','0');

        fclose(fid);
    catch exception
        disp('Ampac_to_xyz failed!');
        throw(exception)
    end

    [~,~] = system(['"', ampacexe, '" "', full_path,'.dat"']);

    outfile = fileread([full_path,'.out']);
    outfile = textscan(outfile,'%s','delimiter','\n');
    outfile = outfile{1};
    
    for i = length(outfile):-1:1
        if (~isempty(regexpi(outfile{i},'Cartesian Coordinates','match')))
            break;
        end
    end
    
    if (i==1)
        throw(MException('AmpacToXYZ:AmpacFailed','Ampac did not run correctly'));
    end
    
    xyz = zeros(natoms,3);
    count = 1;
    
    for j = i+2:i+2+natoms-1
        tmp = regexpi(outfile{j},'\s*\d+\s+\D+\s+(\S+)\s+(\S+)\s+(\S+)\s*$','tokens');
        xyz(count,:) = str2double(tmp{1});
        count = count + 1;
    end
        
    
    atom_type = {zmatrix{:,2}};
    delete([full_path,'.dat']);
    delete([full_path,'.out']);
    delete([full_path,'.vis']);
    delete([full_path,'.arc']);
end

