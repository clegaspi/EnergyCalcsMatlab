function [wl,abs,varargout] = import_uvvis(filename)

fid = fopen(filename,'r');
raw = textscan(fid,'%s','delimiter','\n');
raw = raw{1};

if (nargout > 2)
    varargout{1} = raw{1}(1:end-2);
end

for i = length(raw)-1:-1:3
    if (isempty(raw{i}))
        raw = raw(3:i-1);
        break;
    end
end

wl = zeros(length(raw),1);
abs = wl;

for i = 1:length(raw)
    tmp = textscan(raw{i},'%f %f','delimiter',',');
    wl(i) = tmp{1};
    abs(i) = tmp{2};
end