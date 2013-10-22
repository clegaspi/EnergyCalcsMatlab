function [ rgb ] = color_hex2rgb( hex, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if (length(hex) ~= 6 && length(hex) ~= 7)
        throw(MException('COLOR_HEX2RGB:BadLength','Hex string must have form ''#123456'' or ''123456'''));
    end
    if (isempty(regexpi(hex,'#?[A-F0-9]{6}')) || regexpi(hex,'#?[A-F0-9]{6}') ~= 1)
        throw(MException('COLOR_HEX2RGB:MisformedHex','Hex string must have form ''#123456'' or ''123456'''));
    end
    
    col = regexpi(hex,'#?(?<r>[A-F0-9]{2})(?<g>[A-F0-9]{2})(?<b>[A-F0-9]{2})','names');
    
    rgb = [0 0 0];
    
    rgb(1) = hex2dec(col.r);
    rgb(2) = hex2dec(col.g);
    rgb(3) = hex2dec(col.b);
    
    if (~isempty(varargin))
        if (strcmpi(varargin{1},'normal'))
            return;
        elseif (~strcmpi(varargin{1},'matlab'))
            warning('COLOR_HEX2RGB:UnknownKeyword',['Keyword ''',varargin{1},''' not recognized. Ignoring.']);
        end
    end
    
    rgb = rgb ./ 255;
end

