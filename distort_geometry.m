function newzmatrix = distort_geometry(zmatrix, ring_info, newblsbyatom, varargin)
    AmpacEXE = 'C:\Program Files\Semichem, Inc\Ampac-10.1\ampac.exe';    
    newzmatrix = zmatrix;
    natoms = size(zmatrix, 1);
    rings = {ring_info.rings};
    ringonly = ring_info.ringonly;
    ringinc = ring_info.ringinc;
    
    if (~isempty(varargin) && strcmpi(varargin{1}, 'am1optimize'))
        am1_optimize = true;
    else
        am1_optimize = false;
    end
    
    % Get atom numbers for all angles in z-matrix
    allbas = zeros(natoms,natoms,natoms);
    for i = 3:natoms
        allbas(i, str2double(zmatrix{i,6}), str2double(zmatrix{i,7})) = str2double(zmatrix{i,4});
        allbas(str2double(zmatrix{i,7}), str2double(zmatrix{i,6}), i) = str2double(zmatrix{i,4});
    end

    % Write new bond lengths to new z-matrix
    for i = 2:natoms
        if (strcmpi(zmatrix{i,2},'c'))              
            newzmatrix{i,3} = num2str(newblsbyatom(i,str2double(newzmatrix{i,6})));
        end
    end
    
    if (~am1_optimize)
        newringang = cell(length(rings),1);

        for i = 1:length(newringang)
            blset = zeros(length(rings{i}),1);
            for j = 1:length(rings{i})-1
                blset(j) = newblsbyatom(rings{i}(j), rings{i}(j+1));
            end
            blset(end) = newblsbyatom(rings{i}(end), rings{i}(1));

            ropt = RingOptimizer(blset);
            [in,~,~,~] = ropt.run('deg');
            newringang{i} = in;
        end

        internalringangles = cell(sum(cellfun(@(x)length(x), rings)), 2);

        i = 1;

        for j = 1:length(rings)
            ratoms = [rings{j} rings{j}(1:2)];
            for k = 1:length(ratoms)-2
                internalringangles{i,1} = ratoms(k:(k+2));
                internalringangles{i,2} = newringang{j}(k);
                i = i + 1;
            end
        end

        % Get internal ring angles which match the atoms that connect to
        % rings
        results = cell(length(ringinc),1);
        ira = reshape(cell2mat({internalringangles{:,1}}),3,[])';
        for i = 1:length(ringinc)
            centeratom = ringinc(i,2);

        %     for j = 1:length(allbanums)
        %         if (~all(allbanums(j,:) == ringinc(i,:)))
        %             if (allbanums(j,2) == centeratom)
        %                 results{i}(end+1) = j;
        %             end
        %         end
        %     end

            for j = 1:length(ira)
                if (~all(ira(j,:) == ringinc(i,:)))
                    if (ira(j,2) == centeratom)
                        results{i}(end+1) = j;
                    end
                end
            end
        end

        % Split the internal angle in half and make that the bond angle. If
        % the center atom is part of two rings, subtract the two internal
        % angles from 360 to get the bond angle necessary
        for i = 1:length(results)
            if (length(results{i}) == 1)
                nang = internalringangles{results{i}(1),2};
                nang = 180 - nang / 2;
            else
                nang = internalringangles{results{i}(1),2};
                nang = 360 - nang - internalringangles{results{i}(2),2};
            end
            allbas(ringinc(i,1),ringinc(i,2),ringinc(i,3)) = nang;
            allbas(ringinc(i,3),ringinc(i,2),ringinc(i,1)) = nang;
        end

        % Assign only internal ring angles their new values
        for i = 1:length(ringonly)
            loc = find(cellfun(@(x)isequal(x,ringonly(i,:)), internalringangles(:,1)));
            if (isempty(loc))
                loc = find(cellfun(@(x)isequal(x,fliplr(ringonly(i,:))), internalringangles(:,1)));
            end

            if (isempty(loc))   % This means the angle we want is part of two rings...gotta get the exact value
                % This part of the algorithm assumes the carbons are
                % sequentially numbered...hopefully this won't be an issue
                big_atom_num = max(ringonly(i,[1 3]));
                little_atom_num = min(ringonly(i,[1 3]));

                a1 = find(cellfun(@(x)isequal(x(1:2),[big_atom_num, ringonly(i,2)]), internalringangles(:,1)));
                if (isempty(a1))
                    a1 = find(cellfun(@(x)isequal(x(2:3),[ringonly(i,2), big_atom_num]), internalringangles(:,1)));
                    tmp = internalringangles{a1,1};
                    other_atom = tmp(1);
                else
                    tmp = internalringangles{a1,1};
                    other_atom = tmp(3);
                end

                a2 = find(cellfun(@(x)isequal(x,[little_atom_num, ringonly(i,2), other_atom]), internalringangles(:,1)));
                if (isempty(a2))
                    a2 = find(cellfun(@(x)isequal(x,[other_atom, ringonly(i,2), little_atom_num]), internalringangles(:,1)));
                end

                new_ang = 360 - internalringangles{a1,2} - internalringangles{a2,2};
            else
                new_ang = internalringangles{loc,2};
            end


            allbas(ringonly(i,1),ringonly(i,2),ringonly(i,3)) = new_ang;
            allbas(ringonly(i,3),ringonly(i,2),ringonly(i,1)) = new_ang;
        end

        % Assign new bond angles into the new z-matrix
        for i = 3:natoms
            nang = allbas(i, str2double(newzmatrix{i,6}), str2double(newzmatrix{i,7}));
            newzmatrix{i,4} = num2str(nang);
        end
    else
        try
            params = 'am1 rhf singlet t=auto geo=ok';
            myrand = int32(rand*1e5);
            fid = fopen([pwd,'\distort_geom_temp_',num2str(myrand),'.dat'],'w');

            fprintf(fid, '%s\r\n%s\r\n%s\r\n', params, 'Title', 'Comment');

            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{1,2}, '0.000000', '1', '0.000000', '1', '0.000000', '1', '0','0','0');

            if (natoms > 1)
                fprintf(fid,'%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{2,2}, newzmatrix{2,3}, '0',...
                    '0.000000', '1', '0.000000', '1', newzmatrix{2,6},'0','0');
            end
            if (natoms > 2)
                fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{3,2}, newzmatrix{3,3}, '0',...
                    newzmatrix{3,4}, '1', '0.000000', '1', newzmatrix{3,6}, newzmatrix{3,7},'0');
            end
            if (natoms > 3)
                for i = 4:natoms
                    fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{i,2}, newzmatrix{i,3}, '0',...
                        newzmatrix{i,4}, '1', newzmatrix{i,5}, '1', newzmatrix{i,6}, newzmatrix{i,7}, newzmatrix{i,8});
                end
            end

            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s', '0', '0.000000', '0', '0.000000', '0', '0.000000', '0', '0','0','0');

            fclose(fid);
        catch exception
            disp('Distort geometry failed!');
            throw(exception)
        end
        
        [~,~] = system(['"', AmpacEXE, '" "', pwd, '\distort_geom_temp_',num2str(myrand),'.dat"']);
        
        newzmatrix = ampac_to_zmatrix([pwd,'\distort_geom_temp_',num2str(myrand),'.out']);
        delete([pwd,'\distort_geom_temp_',num2str(myrand),'.dat']);
        delete([pwd,'\distort_geom_temp_',num2str(myrand),'.out']);
        delete([pwd,'\distort_geom_temp_',num2str(myrand),'.vis']);
        delete([pwd,'\distort_geom_temp_',num2str(myrand),'.arc']);
    end
end