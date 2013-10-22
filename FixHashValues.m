path = 'o:\';
% EnergyCalcExp.generate_catalog(path);
% S = load([path, 'catalog.mat']);
% lib = S.lib;
% mykeys = lib.keys;
% myvals = lib.values;

for ikey = 1:length(mykeys)
    curkey = lib(mykeys{ikey});
    efm = curkey.Efield_mag;
    
    newefm = str2double(num2str(efm, 5));
%     newefm = num2str(int32(floor(efm)));
%     efm = efm - str2double(newefm);
%     newefm = [newefm, '.'];
%     
%     decimal_buffer = [];
%     getOut = false;
%     for prec = 1:10
%         if (~getOut)
%             efm = efm * 10;
%             digit = num2str(int32(floor(efm)));
%             if (strcmp(digit,'0') && (length(decimal_buffer) < 5))
%                 decimal_buffer = [decimal_buffer, '0'];
%             elseif (strcmp(digit,'0') && (length(decimal_buffer) == 5))
%                 getOut = true;
%             else
%                 newefm = [newefm, decimal_buffer, digit];
%             end
%             efm = efm - str2double(digit);
%         end
%     end
    
    curkey.Efield_mag = str2double(newefm);
    hash = DataHash(curkey);
    
    if (~strcmp(mykeys{ikey}, hash))
        if (exist([path, mykeys{ikey}, '.mat'], 'file'))
            movefile([path, mykeys{ikey}, '.mat'],[path, hash, '.mat']);
            S = load([path, hash, '.mat']);
            S.key = curkey;
            S.hash = hash;
            save([path, hash, '.mat'], '-struct', 'S');
        end
        
        movefile([path, mykeys{ikey}, '-dm.bin'],[path, hash, '-dm.bin']);
        
        
    end
end
            
        