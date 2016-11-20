function R = project_to_SOd(M)
% function R = project_to_SOd(M)
%
% Given a square d x d matrix M, this function computes and returns a 
% closest matrix R \in SO(d)

d = size(M, 1);

[U, S, V] = svd(M);

R = U * diag([ones(1, d-1), sign(det(U*V'))]) * V';

end

