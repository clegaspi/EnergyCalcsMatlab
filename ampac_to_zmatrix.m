function zmatrix = ampac_to_zmatrix( full_path, varargin )
%AMPAC_TO_ZMATRIX Summary of this function goes here
%   Detailed explanation goes here
    
    if (isempty(varargin) || varargin{1})
        leave_as_strings = true;
    else
        leave_as_strings = false;
    end
    
    if (strcmp(full_path((end-2):end), 'out'))
        outfile = true;
    elseif (strcmp(full_path((end-2):end), 'dat'))
        outfile = false;
    else
        warning('ampac_to_zmatrix:FileExtNotRecognized','ampac_to_zmatrix: Bad file extension');
    end
    
    fdat = fileread(full_path);
    zmat = textscan(fdat,'%s','delimiter','\n');
    zmat = zmat{1};
    
    if (outfile)
        idx = 0;
        i = 1;
        foundonce = false;

        while (idx == 0 && i <= length(zmat))
            %if (strcmp(strtrim(zmat{i}), 'GEOMETRY OPTIMISED : ENERGY MINIMISED'))
                %while (idx == 0)
                    if (~isempty(regexp(strtrim(zmat{i}),'\(I\).+NC$','match')))
                        if (foundonce)
                            idx = i+1;
                            eidx = 0;
                            i = i+2;
                            while (eidx == 0)
                                if (strcmp(zmat{i},''))
                                    eidx = i;
                                else
                                    i = i + 1;
                                end
                            end
                        else
                            foundonce = true;
                            i = i+1;
                        end
                    else
                        i = i+1;
                    end
                %end
            %else
                %i = i + 1;
            %end
        end
    else
        idx = 4;
        if (isempty(zmat{end}))
            eidx = size(zmat, 1) - 1;
        else
            eidx = size(zmat, 1);
        end
    end
    
    natoms = eidx - idx;
    zmat = {zmat{idx:eidx-1}};
    zmatrix = cell(natoms, 8);
    
    if (outfile)
        temp = textscan(zmat{1}, '%s');
        temp = temp{1};
        zmatrix{1,1} = temp{1};
        zmatrix{1,2} = temp{2};
        zmatrix{1,3} = '';
        zmatrix{1,4} = '';
        zmatrix{1,5} = '';
        zmatrix{1,6} = '';
        zmatrix{1,7} = '';
        zmatrix{1,8} = '';

        temp = textscan(zmat{2}, '%s');
        temp = temp{1};
        zmatrix{2,1} = temp{1};
        zmatrix{2,2} = temp{2};
        zmatrix{2,3} = temp{3};
        zmatrix{2,4} = '';
        zmatrix{2,5} = '';
        
        if (strcmp(temp{4},'*'))
            zmatrix{2,6} = temp{5};
        else
            zmatrix{2,6} = temp{4};
        end
            
        zmatrix{2,7} = '';
        zmatrix{2,8} = '';

        temp = textscan(zmat{3}, '%s');
        temp = temp{1};
        zmatrix{3,1} = temp{1};
        zmatrix{3,2} = temp{2};
        zmatrix{3,3} = temp{3};
        if (strcmp(temp{4},'*'))
            lastidx = 5;
            zmatrix{3,4} = temp{5};
        else
            lastidx = 4;
            zmatrix{3,4} = temp{4};
        end
        if (strcmp(temp{lastidx+1},'*'))
            lastidx = lastidx+2;
        else
            lastidx = lastidx+1;
        end
        zmatrix{3,5} = '';
        zmatrix{3,6} = temp{lastidx};
        zmatrix{3,7} = temp{lastidx+1};
        zmatrix{3,8} = '';

        for i = 4:natoms
            temp = textscan(zmat{i}, '%s');
            temp = temp{1};
            zmatrix{i,1} = temp{1};
            zmatrix{i,2} = temp{2};
            zmatrix{i,3} = temp{3};
            if (strcmp(temp{4},'*'))
                lastidx = 5;
                zmatrix{i,4} = temp{5};
            else
                lastidx = 4;
                zmatrix{i,4} = temp{4};
            end
            if (strcmp(temp{lastidx+1},'*'))
                zmatrix{i,5} = temp{lastidx+2};
                lastidx = lastidx+2;
            else
                zmatrix{i,5} = temp{lastidx+1};
                lastidx = lastidx+1;
            end
            if (strcmp(temp{lastidx+1},'*'))
                zmatrix{i,6} = temp{lastidx+2};
                lastidx = lastidx+2;
            else
                zmatrix{i,6} = temp{lastidx+1};
                lastidx = lastidx+1;
            end
            zmatrix{i,7} = temp{lastidx+1};
            zmatrix{i,8} = temp{lastidx+2};
        end
    else
        temp = textscan(zmat{1}, '%s');
        temp = temp{1};
        zmatrix{1,1} = '1';
        zmatrix{1,2} = temp{1};
        zmatrix{1,3} = '';
        zmatrix{1,4} = '';
        zmatrix{1,5} = '';
        zmatrix{1,6} = '';
        zmatrix{1,7} = '';
        zmatrix{1,8} = '';

        temp = textscan(zmat{2}, '%s');
        temp = temp{1};
        zmatrix{2,1} = '2';
        zmatrix{2,2} = temp{1};
        zmatrix{2,3} = temp{2};
        zmatrix{2,4} = '';
        zmatrix{2,5} = '';
        zmatrix{2,6} = temp{8};
        zmatrix{2,7} = '';
        zmatrix{2,8} = '';

        temp = textscan(zmat{3}, '%s');
        temp = temp{1};
        zmatrix{3,1} = '3';
        zmatrix{3,2} = temp{1};
        zmatrix{3,3} = temp{2};
        zmatrix{3,4} = temp{4};
        zmatrix{3,5} = '';
        zmatrix{3,6} = temp{8};
        zmatrix{3,7} = temp{9};
        zmatrix{3,8} = '';

        for i = 4:natoms
            temp = textscan(zmat{i}, '%s');
            temp = temp{1};
            zmatrix{i,1} = num2str(i);
            zmatrix{i,2} = temp{1};
            zmatrix{i,3} = temp{2};
            zmatrix{i,4} = temp{4};
            zmatrix{i,5} = temp{6};
            zmatrix{i,6} = temp{8};
            zmatrix{i,7} = temp{9};
            zmatrix{i,8} = temp{10};
        end
    end
    
    if (~leave_as_strings)
        for i = 1:natoms
            for j = [1,3:8]
                if (~isempty(zmatrix{i,j}))
                    zmatrix{i,j} = str2double(zmatrix{i,j});
                end
            end
        end
    end

end

