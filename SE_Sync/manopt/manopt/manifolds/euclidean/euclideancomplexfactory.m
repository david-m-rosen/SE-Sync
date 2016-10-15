function M = euclideancomplexfactory(m, n)
% Returns a manifold struct to optimize over complex m-by-n matrices.
%
% function M = euclideancomplexfactory(m, n)
%
% Returns M, a structure describing the vector space of complex m-by-n
% matrices, as a manifold for Manopt.
%
% The complex plane is here viewed as R^2. The inner product between two
% m-by-n matrices A and B is given by: real(trace(A'*B)). This choice
% guides the proper definition of gradient and Hessian for this geometry.
% This is not the classical Euclidean inner product for complex matrices;
% it is a real inner product.
%
% See also: euclideanfactory

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, April 7, 2015.
% Contributors: 
% Change log: 

    
    if ~exist('n', 'var') || isempty(n)
        n = 1;
    end

    M.name = @() sprintf('Vector space C^(%dx%d)', m, n);
    
    M.dim = @() 2*m*n;
    
    M.inner = @(x, d1, d2) real(d1(:)'*d2(:));
    
    M.norm = @(x, d) norm(d, 'fro');
    
    M.dist = @(x, y) norm(x-y, 'fro');
    
    M.typicaldist = @() sqrt(m*n);
    
    M.proj = @(x, d) d;
    
    M.egrad2rgrad = @(x, g) g;
    
    M.ehess2rhess = @(x, eg, eh, d) eh;
    
    M.tangent = M.proj;
    
    M.exp = @exp;
    function y = exp(x, d, t)
        if nargin == 3
            y = x + t*d;
        else
            y = x + d;
        end
    end
    
    M.retr = M.exp;
	
	M.log = @(x, y) y-x;

    M.hash = @(x) ['z' hashmd5([real(x(:)) ; imag(x(:))])];
    
    M.rand = @() (randn(m, n) + 1i*randn(m, n))/sqrt(2);
    
    M.randvec = @randvec;
    function u = randvec(x) %#ok<INUSD>
        u = randn(m, n) + 1i*randn(m, n);
        u = u / norm(u, 'fro');
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(x) zeros(m, n);
    
    M.transp = @(x1, x2, d) d;
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    mn = m*n;
    M.vec = @(x, u_mat) [real(u_mat(:)) ; imag(u_mat(:))];
    M.mat = @(x, u_vec) reshape(u_vec(1:mn), [m, n]) + 1i*reshape(u_vec((mn+1):end), [m, n]);
    M.vecmatareisometries = @() true;

end
