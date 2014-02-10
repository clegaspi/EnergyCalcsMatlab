function [ varargout ] = closest_member( a, myset )
%CLOSEST_MEMBER Finds the member of myset which is closest in value to a
%   INPUT
%   a - scalar, search value
%   myset - vector
%   
%   OUTPUT
%   varargout{1} - closest value of myset to a. If a is equidistant between
%       two values within myset, it will be the left-hand value as sorted
%       ascending (left<right)
%   varargout{2} - right-hand value of myset closest to a. If the
%       varargout{1}==a exactly, then varargout{2}==varargout{1}==a
%   varargout{3} - indices of left- and right-hand values

    [b,ix] = sort(myset,'ascend');
    
    for i = 1:length(b)
        if (b(i) < a && i < length(b))
            continue;
        elseif (b(i) == a || (b(i) < a && i == length(b)))
            varargout{1} = b(i);
            switch nargout
                case 2
                    varargout{2} = b(i);
                case 3
                    varargout{2} = b(i);
                    varargout{3} = [ix(i) ix(i)];
            end
            break;
        else
            varargout{1} = b(i-1);
            switch nargout
                case 2
                    varargout{2} = b(i);
                case 3
                    varargout{2} = b(i);
                    varargout{3} = [ix(i-1) ix(i)];
            end
            break;
        end
    end
end

