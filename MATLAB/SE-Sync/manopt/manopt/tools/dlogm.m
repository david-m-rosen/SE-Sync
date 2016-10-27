function D = dlogm(X, H)
% Fréchet derivative of the matrix logarithm.
%
% function D = dlogm(X, H)
%
% Computes the directional derivative (the Fréchet derivative) of logm at X
% along H (square matrices).
%
% Thus, D = lim_(t -> 0) (logm(X + tH) - logm(X)) / t.
% 
% See also: dfunm dexpm dsqrtm

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2015.
% Contributors:
% Change log:
    
    D = dfunm(@logm, X, H);
    
end
