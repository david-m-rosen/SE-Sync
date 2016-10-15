function M = stiefelgeneralizedfactory(n, p, B)
% Returns a manifold structure of "scaled" orthonormal matrices.
%
% function M = stiefelgeneralizedfactory(n, p)
% function M = stiefelgeneralizedfactory(n, p, B)
%
% The generalized Stiefel manifold is the set of "scaled" orthonormal 
% nxp matrices X such that X'*B*X is identity. B must be positive definite.
% If B is identity, then this is the standard Stiefel manifold.
%
% The generalized Stiefel manifold is endowed with a scaled metric
% by making it a Riemannian submanifold of the Euclidean space,
% again endowed with the scaled inner product.
%
% Some notions (not all) are from Section 4.5 of the paper
% "The geometry of algorithms with orthogonality constraints",
% A. Edelman, T. A. Arias, S. T. Smith, SIMAX, 1998.
%
% Paper link: http://arxiv.org/abs/physics/9806030.
%
% Note: egrad2rgrad and ehess2rhess involve solving linear systems in B. If
% this is a bottleneck for a specific application, then a way forward is to
% create a modified version of this file which preprocesses B to speed this
% up (typically, by computing a Cholesky factorization of it, then calling
% an appropriate solver).
%
% See also: stiefelfactory  grassmannfactory  grassmanngeneralizedfactory 

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, June 30, 2015.
% Contributors:
%
% Change log:
%   

    
    if ~exist('B', 'var') || isempty(B)
        B = speye(n); % Standard Stiefel manifold.
    end
    
    M.name = @() sprintf('Generalized Stiefel manifold St(%d, %d)', n, p);
    
    M.dim = @() (n*p - .5*p*(p+1));
    
    M.inner = @(X, eta, zeta) trace(eta'*(B*zeta)); % Scaled metric.
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(X, Y) error('stiefelgeneralizedfactory.dist not implemented yet.');
    
    M.typicaldist = @() sqrt(p);
    
    % Orthogonal projection of an ambient vector U to the tangent space
    % at X.
    M.proj = @projection;
    function Up = projection(X, U)
        BX = B*X;
        
        % Projection onto the tangent space
        Up = U - X*symm(BX'*U);  
    end
    
    M.tangent = M.proj;
    
    M.egrad2rgrad = @egrad2rgrad;
    function rgrad = egrad2rgrad(X, egrad)
        
        % First, scale egrad according the to the scaled metric in the
        % Euclidean space.
        egrad_scaled = B\egrad;
        
        % Second, project onto the tangent space.
        % rgrad = egrad_scaled - X*symm((B*X)'*egrad_scaled);
        %
        % Verify that symm(BX'*egrad_scaled) = symm(X'*egrad).
        
        rgrad = egrad_scaled - X*symm(X'*egrad);
    end
    
    
    
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, H)
        egraddot = ehess;
        Xdot = H;
        
        % Directional derivative of the Riemannian gradient.
        egrad_scaleddot = B\egraddot;
        rgraddot = egrad_scaleddot - Xdot*symm(X'*egrad)...
            - X*symm(Xdot'*egrad)...
            - X*symm(X'*egraddot);
        
        % Project onto the tangent space.
        rhess = M.proj(X, rgraddot);
    end
    
    
    M.retr = @retraction;
    function Y = retraction(X, U, t)
        if nargin < 3
            t = 1.0;
        end
        Y = guf(X + t*U); % Ensure that Y'*B*Y is identity.
    end
    
    
    M.exp = @exponential;
    function Y = exponential(X, Z, t)
        if nargin < 3
            t = 1.0;
        end
        Y = retraction(X, Z, t);
        warning('manopt:stiefelgeneralizedfactory:exp', ...
               ['Exponential for generalized Stiefel manifold ' ...
                'manifold not implemented yet. Used retraction instead.']);
    end


    M.hash = @(X) ['z' hashmd5(X(:))];
    
    M.rand = @random;
    function X = random()
        X = guf(randn(n, p)); % Ensure that X'*B*X is identity;
    end
    
    M.randvec = @randomvec;
    function U = randomvec(X)
        U = projection(X, randn(n, p));
        U = U / norm(U(:));
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(X) zeros(n, p);
    
    % This transport is compatible with the generalized polar retraction.
    M.transp = @(X1, X2, d) projection(X2, d);
    
    M.vec = @(X, u_mat) u_mat(:);
    M.mat = @(X, u_vec) reshape(u_vec, [n, p]);
    M.vecmatareisometries = @() false;
    
    % Some auxiliary functions
    symm = @(D) (D + D')/2;
    
    function X = guf(Y)
        % Generalized polar decomposition of an n-by-p matrix Y.
        % X'*B*X is identity.
        
        % Method 1
        [u, ~, v] = svd(Y, 0);
  
        % Instead of the following three steps, an equivalent, but an 
        % expensive way is to do X = u*(sqrtm(u'*(B*u))\(v')).
        [q, ssquare] = eig(u'*(B*u));
        qsinv = q/sparse(diag(sqrt(diag(ssquare))));
        X = u*((qsinv*q')*v'); % X'*B*X is identity.
        
        
        % Another computation using restricted_svd
        % [u, ~, v] = restricted_svd(Y);
        % X = u*v'; % X'*B*X is identity.
        
    end
    
    function [u, s, v] = restricted_svd(Y)
        % We compute a thin svd-like decomposition of an n-by-p matrix Y 
        % into matrices u, s, and v such that u is an n-by-p matrix
        % with u'*B*u being identity, s is a p-by-p diagonal matrix 
        % with positive entries, and v is a p-by-p orthogonal matrix.
        % Y = u*s*v'.
        [v, ssquare] = eig(symm(Y'*(B*Y))); % Y*B*Y is positive definite
        ssquarevec = diag(ssquare);
        
        s = sparse(diag(abs(sqrt(ssquarevec))));
        u = Y*(v/s); % u'*B*u is identity.
    end

end
