function N = tangentspherefactory(M, x)
% Returns a manifold struct. for the sphere on the tangent space to M at x.
%
% N = tangentspherefactory(M, x)
%
% N defines a manifold that is the unit sphere on the tangent space to M
% at x. Points are represented as tangent vectors of unit norm. Tangent
% vectors are represented as tangent vectors orthogonal to the root point,
% with respect to the Riemannian metric on the tangent space.
%
% This is chiefly useful to solve optimization problems involving unit norm
% tangent vectors to M at x, which notably comes up when looking for
% extreme eigenvectors of the Hessian of a cost function on M at x, for
% example. The Riemannian structure on this sphere is that of a Riemannian
% submanifold of the (Euclidean) tangent space, equipped with the
% Riemannian metric of M at that point.
%
% See also: hessianextreme

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, March 16, 2015.
% Contributors: 
% Change log: 

    % N is the manifold we build. y will be a point on N, thus also a
    % tangent vector to M at x. This is a typical Riemannian submanifold of
    % a Euclidean space, hence it will be easy to describe in terms of the
    % tools available for M.
    N = struct();
    
    % u, u1 and u2 will be tangent vectors to N at y. The tangent space to
    % N at y is a subspace of the tangent space to M at x, thus u, u1 and
    % u2 are also tangent vectors to M at x.
    
    N.dim   = @() M.dim() - 1;
    N.inner = @(y, u1, u2) M.inner(x, u1, u2);
    N.norm  = @(y, u)      M.norm(x, u);
    N.proj  = @(y, v) M.lincomb(x, 1, v, -M.inner(x, v, y), y);
    N.typicaldist = @() 1;
    N.tangent = N.proj;
    N.egrad2rgrad = N.proj;
    N.retr = @retraction;
    N.exp = N.retr;
    function yy = retraction(y, u, t)
        if nargin == 2
            t = 1;
        end
        y_plus_tu = M.lincomb(x, 1, y, t, u);
        nrm = M.norm(x, y_plus_tu);
        yy = M.lincomb(x, 1/nrm, y_plus_tu);
    end
    N.rand = @random;
    function y = random()
        y = M.randvec(x);
        nrm = M.norm(x, y);
        y = M.lincomb(x, 1/nrm, y);
    end
    N.randvec = @randvec;
    function u = randvec(y)
        u = N.proj(y, N.rand());
        nrm = N.norm(y, u);
        u = M.lincomb(x, 1/nrm, u);
    end
    N.zerovec = M.zerovec;
    N.lincomb = M.lincomb;
    N.transp = @(y1, y2, u) N.proj(y2, u);
    N.hash = @(y) ['z' hashmd5(M.vec(x, y))];
    
end
