function [ hex ] = color_rgb2hex( rgb, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if (numel(rgb) ~= 3)
        throw(MException('COLOR_RGB2HEX:BadLength','RGB vector must be 1x3 or 3x1'));
    end
    matlab_type = false;   % Normal where max is 255
    if (isempty(varargin) || strcmpi(varargin{1},'matlab'))
        rgb = rgb .* 255;
        matlab_type = true;   % Matlab RGB where max is 1
    end
    if (any(rgb > 255))
        max_val = '255';
        if (matlab_type)
            max_val = '1';
        end
        throw(MException('COLOR_RGB2HEX:MisformedHex',['All RGB values must be <= ',max_val,' for specified RGB type']));
    end
    rgb = int32(rgb);
    hex = [dec2hex(rgb(1),2), dec2hex(rgb(2),2), dec2hex(rgb(3),2)];
end

