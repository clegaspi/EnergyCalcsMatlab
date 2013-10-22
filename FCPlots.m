% To read in data, export from origin/excel as a csv file. Drag file into
% matlab command window and go through prompts. On the second screen, you
% have the option to change the name of the variable. Change it to "data"
% so that the script will continue to function. If you forget to do that,
% you can take whatever variable name it gave it (say it called it "test"
% and run the following line in the command window

% data = test;

% Without the percent sign. The percent sign is the comment character so
% matlab will ignore all lines of code with a percent sign. The double
% percent is the marker of the beginning of a new cell.

% Remember to only convert the data once! Otherwise, you will be
% flip-flopping your energy/wavelength values and your graph will look
% weird!

% To label your axes and title the graph, type the following into the
% command window without the % signs:

% title('Insert title of graph here');
% xlabel('Insert x label here');
% ylabel('Insert y label here');

% If you have more than one figure window open, make sure you specify which
% figure window to plot in by either clicking on the plot window or typing
% figure(n) where n is the figure number into the command window.
% You can make numbers superscripted by inserting a ^ before the character
% in the label.

%% Read in data and convert
data(:,1) = 1e7 ./ data(:,1);   % changing to wavenums
data(:,2) = data(:,2) / max(data(:,2)); % normalizing y values

%% 
set(gca,'xdir','reverse');  % Reverses x axis to be descending
figure(1);
hold on;

plot(data(:,1),data(:,2));  % plotting

[maxdat,~] = ginput(1); % selecting 0-0 peak

%% Generate FC Plot
delta = 0.8;
qn = 0:10;
CCstretch = 1300;

intensity = ((0.5 * delta ^ 2) .^ qn) ./ factorial(qn); % Calculate intensity
figure(1);
hold on;    % Keeps previous graph on the plot
bar(maxdat - qn * CCstretch, intensity, 'BarWidth', 0.25);  % Plots bars as function of energy