function D = dsqrtm(X, H)
% Fréchet derivative of the matrix square root.
%
% function D = dsqrtm(X, H)
%
% Computes the directional derivative (the Fréchet derivative) of sqrtm at
% X along H (square matrices).
%
% Thus, D = lim_(t -> 0) (sqrtm(X + tH) - sqrtm(X)) / t.
% 
% See also: dfunm dlogm dexpm

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2015.
% Contributors:
% Change log:
    
    D = dfunm(@sqrtm, X, H);
    
end
