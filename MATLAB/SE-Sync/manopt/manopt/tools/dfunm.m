function D = dfunm(funm, X, H)
% Fréchet derivative of matrix functions.
%
% function D = dfunm(funm, X, H)
%
% Computes the directional derivative (the Fréchet derivative) of a matrix
% function (such as @logm, @expm, ...) at X along H (square matrices),
% according to a very nice trick which appears in this paper:
% 
% "Computing the Fréchet derivative of the matrix exponential, with an
% application to condition number estimation",
% Awad H. Al-Mohy and Nicholas J. Higham, 2009.
% http://eprints.ma.man.ac.uk/1218/01/covered/MIMS_ep2008_26.pdf
%
% Thus, D = lim_(t -> 0) (funm(X + tH) - funm(X)) / t.
%
% This code is simple, but may not be the most efficient. In particular, it
% requires computing the matrix function on matrices which are four times
% as big, and which may have lost important structure (such as symmetry).
% 
% See also: dlogm dexpm dsqrtm

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2015.
% Contributors:
% Change log:
    
    n = size(X, 1);
    
    assert(length(size(X)) == 2,     'X and H must be square matrices.');
    assert(length(size(H)) == 2,     'X and H must be square matrices.');
    assert(size(X, 1) == size(X, 2), 'X and H must be square matrices.');
    assert(all(size(X) == size(H)),  'X and H must have the same size.');
    
    Z = zeros(n);
    A = funm([X, H ; Z, X]);
    D = A(1:n, (n+1):end);
    
end
