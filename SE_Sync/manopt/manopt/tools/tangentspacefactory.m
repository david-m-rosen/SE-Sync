function N = tangentspacefactory(M, x)
% Returns a manifold structure representing the tangent space to M at x.
%
% N = tangentspacefactory(M, x)
%
% N defines a (linear) manifold that is the tangent space to M at x. Points
% are represented as tangent vectors to M at x. Tangent vectors are also
% represented as tangent vectors to M at x.
%
% This is chiefly useful to solve optimization problems involving tangent
% vectors to M at x, which notably comes up when solving linear systems
% involving, for example, the Hessian of the cost on M at x. The Riemannian
% (actually, Euclidean) structure on N is that of the tangent space to M,
% that is, the inner product is inherited.
%
% See also: preconhessiansolve

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, April 9, 2015.
% Contributors: 
% Change log: 

    % N is the manifold we build. y will be a point on N, thus also a
    % tangent vector to M at x. This is a typical Euclidean space, hence it
    % will be easy to describe in terms of the tools available for M.
    N = struct();
    
    % u, u1 and u2 will be tangent vectors to N at y. The tangent space to
    % N at y is the tangent space to M at x, thus u, u1 and u2 are also
    % tangent vectors to M at x.
    
    N.dim   = @() M.dim();
    N.inner = @(y, u1, u2) M.inner(x, u1, u2);
    N.norm  = @(y, u) M.norm(x, u);
    N.proj  = M.proj;
    N.typicaldist = @() N.dim();
    N.tangent = @(y, u) u;
    N.egrad2rgrad = @(x, g) g;
    N.ehess2rhess = @(x, eg, eh, d) eh;
    N.exp = @exponential;
    N.retr = @exponential;
    N.log = @(y1, y2) M.lincomb(x, 1, y2, -1, y1);
    N.pairmean = @(y1, y2) M.lincomb(x, 0.5, y1, 0.5, y2);
    N.rand = @() M.randvec(x);
    N.randvec = @(y) M.randvec(x);
    N.zerovec = M.zerovec;
    N.lincomb = M.lincomb;
    N.transp = @(y1, y2, u) u;
    N.hash = @(y) ['z' hashmd5(M.vec(x, y))];
    
    % In a Euclidean space, the exponential is merely the sum: y + tu.
    function yy = exponential(y, u, t)
        if nargin == 2
            t = 1;
        end
        yy = M.lincomb(x, 1, y, t, u);
    end
    
end
