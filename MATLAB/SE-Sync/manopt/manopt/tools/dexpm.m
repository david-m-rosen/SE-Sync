function D = dexpm(X, H)
% Fréchet derivative of the matrix exponential.
%
% function D = dexpm(X, H)
%
% Computes the directional derivative (the Fréchet derivative) of expm at X
% along H (square matrices).
%
% Thus, D = lim_(t -> 0) (expm(X + tH) - expm(X)) / t.
% 
% See also: dfunm dlogm dsqrtm

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2015.
% Contributors:
% Change log:
    
    D = dfunm(@expm, X, H);
    
end
