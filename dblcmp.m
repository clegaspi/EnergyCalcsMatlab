function [ res ] = dblcmp( left, right, varargin )
%DBLCMP Summary of this function goes here
%   Detailed explanation goes here
    
    res = false;
    
    if (int32(log10(left)) ~= int32(log10(right)))
        return;
    end
    
    tolerance = left / 1e15;
    
    if (tolerance == 0)
        tolerance = right / 1e15;
        if (tolerance == 0)
            tolerance = 1e-15;
        end
    end
    
    if (~isempty(varargin))
        tolerance = varargin{1};
    end
    
    if (left > right - tolerance && left < right + tolerance)
        res = true;
    end
    
end