function X = removeNaN(X, substitute)

% REMOVENAN Remove NaNs.
% FORMAT
% DESC Substitute NaN entries with the given character argument.
% ARG X : The vector to search for NaNs.
% ARG substitute : The characters to substitute for NaNs.
% RETURN X : Vector with NaNs removed.
% 
% COPYRIGHT : Alfredo Kalaitzis, 2010, 2011
%
% GPREGE

if nargin == 1
    substitute = 0;
end

X(isnan(X)) = substitute;
