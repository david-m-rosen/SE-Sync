function M = multinomialfactory(n, m)
% Manifold of n-by-m column-stochastic matrices with positive entries.
%
% function M = multinomialfactory(n, m)
%
% The returned structure M is a Manopt manifold structure to optimize over
% the set of n-by-m matrices with (strictly) positive entries and such that
% the entries of each column sum to one.
%
% The metric imposed on the manifold is the Fisher metric such that 
% the set of n-by-m column-stochastic matrices (aka the multinomial manifold)
% is a Riemannian submanifold of the space of n-by-m matrices. Also it
% should be noted that the retraction operation that we define 
% is first order and as such the checkhessian tool cannot verify 
% the slope correctly.
%             
% The file is based on developments in the research paper
% Y. Sun, J. Gao, X. Hong, B. Mishra, and B. Yin,
% "Heterogeneous tensor decomposition for clustering via manifold
% optimization", arXiv:1504.01777, 2015.
%
% Link to the paper: http://arxiv.org/abs/1504.01777.
%
% Please cite the Manopt paper as well as the research paper:
%     @Techreport{sun2014multinomial,
%       Title   = {Heterogeneous tensor decomposition for clustering via manifold optimization},
%       Author  = {Sun, Y. and Gao, J. and Hong, X. and Mishra, B. and Yin, B.},
%       Journal = {Arxiv preprint arXiv:1504.01777},
%       Year    = {2014}
%     }

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, April 06, 2015.
% Contributors:
% Change log:
    
    M.name = @() sprintf('%dx%d column-stochastic matrices with positive entries', n, m);
    
    M.dim = @() (n-1)*m;
    
    % We impose the Fisher metric.
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        ip = sum((eta(:).*zeta(:))./X(:));
    end
    
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(X, Y) error('multinomialfactory.dist not implemented yet.');
    
    M.typicaldist = @() m*pi/2; % This is an approximation.
    
    % Column vector of ones of length n. 
    e = ones(n, 1);
    
    M.egrad2rgrad = @egrad2rgrad;
    function rgrad = egrad2rgrad(X, egrad)
        lambda = -sum(X.*egrad, 1); % Row vector of length m.
        rgrad = X.*egrad + (e*lambda).*X; % This is in the tangent space.
    end
    
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, eta)
        
        % Riemannian gradient computation.
        % lambda is a row vector of length m.
        lambda = - sum(X.*egrad, 1);
        rgrad =  X.*egrad + (e*lambda).*X;
        
        % Directional derivative of the Riemannian gradient.
        % lambdadot is a row vector of length m.
        lambdadot = -sum(eta.*egrad, 1) - sum(X.*ehess, 1); 
        rgraddot = eta.*egrad + X.*ehess + (e*lambdadot).*X + (e*lambda).*eta;
        
        % Correction term because of the non-constant metric that we
        % impose. The computation of the correction term follows the use of
        % Koszul formula.
        correction_term = - 0.5*(eta.*rgrad)./X;
        rhess = rgraddot + correction_term;
        
        % Finally, projection onto the tangent space.
        rhess = M.proj(X, rhess);
    end
    
    % Projection of the vector eta in the ambeint space onto the tangent
    % space.
    M.proj = @projection;
    function etaproj = projection(X, eta)
        alpha = sum(eta, 1); % Row vector of length m.
        etaproj = eta - (e*alpha).*X;
    end
    
    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        % A first-order retraction.
        Y = X.*exp(t*(eta./X)); % Based on mapping for positive scalars.
        Y = Y./(e*(sum(Y, 1))); % Projection onto the constraint set.
        % For numerical reasons, so that we avoid entries going to zero:
        Y = max(Y, eps);
    end
    
    M.exp = @exponential;
    function Y = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y = retraction(X, eta, t);
        warning('manopt:multinomialfactory:exp', ...
            ['Exponential for the Multinomial manifold' ...
            'manifold not implemented yet. Used retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5(X(:))];
    
    M.rand = @random;
    function X = random()
        % A random point in the ambient space.
        X = rand(n, m); %
        X = X./(e*(sum(X, 1)));
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(X)
        % A random vector in the tangent space
        eta = randn(n, m);
        eta = M.proj(X, eta); % Projection onto the tangent space.
        nrm = M.norm(X, eta);
        eta = eta / nrm;
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(X) zeros(n, m);
    
    M.transp = @(X1, X2, d) projection(X2, d);
    
    % vec and mat are not isometries, because of the scaled metric.
    M.vec = @(X, U) U(:);
    M.mat = @(X, u) reshape(u, n, m);
    M.vecmatareisometries = @() false;
end
