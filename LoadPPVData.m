% Load PPV data

%% 6-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\6-merPPV\';
mat_file_path = {'..\6-merPPV-EvsF.mat'};

%% 7-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\7-merPPV\';
mat_file_path = {'..\7-merPPV-EvsF.mat'};

%% 8-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\8-merPPV\';
mat_file_path = {'..\8-merPPV-EvsF.mat'};

%% 9-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\9-merPPV\';
mat_file_path = {'..\9-merPPV-EvsF.mat'};

%% 10-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\10-merPPV\';
mat_file_path = {'..\10-merPPV-EvsF.mat'};

%% 13-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\13-merPPV\';
mat_file_path = {'..\13-merPPV-EvsF.mat'};

%% 15-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\15-merPPV\';
mat_file_path = {'..\15-merPPV-EvsF.mat'};


%% 17-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\17-merPPV\';
mat_file_path = {'..\17-merPPV-EvsF.mat'};

%% 20-mer PPV Unrelaxed

loadrootdir = 'C:\Users\clegaspi\Documents\MATLAB\data\20-merPPV\';
mat_file_path = {'..\20-merPPV-EvsF.mat'};

%% Load data here

if (loadrootdir(end) ~= '\')
    loadrootdir = [loadrootdir, '\'];
end

for i = 1:length(mat_file_path)
    load([loadrootdir, mat_file_path{i}]);
end

if (strcmpi(myexp.ampacdatapath, 'p:\'))
    expt = 's';
    indo = 'o';
    ampac = 'p';
else
    expt = 'j';
    indo = 'k';
    ampac = 'l';
end

SFolder = 'Exp';
OFolder = 'INDOLib';
PFolder = 'GSLib';

system(['subst ',expt,': /d']);
system(['subst ',indo,': /d']);
system(['subst ',ampac,': /d']);

system(['subst ',expt,': "',loadrootdir,SFolder,'"']);
system(['subst ',indo,': "',loadrootdir,OFolder,'"']);
system(['subst ',ampac,': "',loadrootdir,PFolder,'"']);