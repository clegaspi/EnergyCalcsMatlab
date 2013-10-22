mols = {'8-merPPV',0.2,[6 62]};
phi1 = [90:5:175];

for imol = 1:size(mols, 1);
    system('subst s: /d');
    system('subst o: /d');
    system('subst p: /d');
    
    rootdir = ['C:\Documents and Settings\Linda Group\My Documents\MATLAB\polymer_calcs\data\',mols{imol,1},'\'];
    SFolder = 'Exp';
    OFolder = 'INDOLib';
    PFolder = 'GSLib';

    system(['subst s: "',rootdir,SFolder,'"']);
    system(['subst o: "',rootdir,OFolder,'"']);
    system(['subst p: "',rootdir,PFolder,'"']);
    
    field = 0:0.001:mols{imol,2};
        
    for k = 1:numel(field)
    %     waithdldata = get(waithdl,'userdata');
    %     if (waithdldata.bail == 0)
    %         waitbar((k-1 + (j-1)*numel(field) + (i-1)*numel(p2i)*numel(field)) / (numel(p1i)*numel(p2i)*numel(field)), waithdl);
    %     else
    %         keyboard;
    %     end
    % 
    %     if (exist(['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat'], 'file'))
    %         S = load(['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat']);
    %         myexp = S.obj;
    %     else
            myece = ECEParams('AM1', {phi1}, true, field(k), 800, 25,...
                's:\1Angle\params.txt');

            if (field(k) ~= 0)
                myece.dm_guess = ['o:\', myexp.data(1).indo_hash, '-dm.bin'];
                myece.try_default_first = true;
            end

            myexp = EnergyCalcExp(myece,...
                ['s:\1Angle\',mols{imol,1},'.dat'],...
                mols{imol,3},...
                ['s:\1Angle\',mols{imol,1},'-',num2str(field(k)),'VA.mat'],...
                'o:\',...
                'p:\',...
                false,...
                ['Bond Order calcs - ',num2str(field(k)),'VA']);
            myexp.run('quiet');
    %     end
    end
end