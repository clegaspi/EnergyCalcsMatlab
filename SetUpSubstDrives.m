%% Unload substituted drives

system('subst s: /d');
system('subst o: /d');
system('subst p: /d');

% If "Invalid parameter" error, then the drives didn't exist as subst drives.

%% Load substituted drives

rootdir = fullfile(pwd,'data\8-merPPV-newampac\Coupling');
SFolder = 'Exp';
OFolder = 'INDOLib';
PFolder = 'GSLib';

if (rootdir(end) ~= '\')
    rootdir = [rootdir, '\'];
end

if (~exist(rootdir,'dir'))
    mkdir(rootdir);
end
if (~exist([rootdir,SFolder],'dir'))
    mkdir([rootdir,SFolder]);
end
if (~exist([rootdir,OFolder],'dir'))
    mkdir([rootdir,OFolder]);
end
if (~exist([rootdir,PFolder],'dir'))
    mkdir([rootdir,PFolder]);
end

system(['subst s: "',rootdir,SFolder,'"']);
system(['subst o: "',rootdir,OFolder,'"']);
system(['subst p: "',rootdir,PFolder,'"']);

%% Unload substituted drives

% SECOND SET IN CASE NEEDED

system('subst j: /d');
system('subst k: /d');
system('subst l: /d');

% If "Invalid parameter" error, then the drives didn't exist as subst drives.

%% Load substituted drives

rootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\MeLPPP-13mer\';
SFolder = 'Exp';
OFolder = 'INDOLib';
PFolder = 'GSLib';

if (rootdir(end) ~= '\')
    rootdir = [rootdir, '\'];
end

if (~exist(rootdir,'dir'))
    mkdir(rootdir);
end
if (~exist([rootdir,SFolder],'dir'))
    mkdir([rootdir,SFolder]);
end
if (~exist([rootdir,OFolder],'dir'))
    mkdir([rootdir,OFolder]);
end
if (~exist([rootdir,PFolder],'dir'))
    mkdir([rootdir,PFolder]);
end

system(['subst j: "',rootdir,SFolder,'"']);
system(['subst k: "',rootdir,OFolder,'"']);
system(['subst l: "',rootdir,PFolder,'"']);