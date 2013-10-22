format long;
cd = fileread('c:\20-merPPV-out.txt');
cd = textscan(cd,'%s','delimiter','\n');
cd = cd{1};
cd = cellfun(@(x)textscan(x,'%s'), cd);

diff = cellfun(@(x)str2double(x{11}), cd);
diff = log10(diff);

figure(54);
plot(1:length(diff), diff);