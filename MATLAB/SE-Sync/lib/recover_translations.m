function t = recover_translations(R, problem_data)
% function t = recover_translations(R, problem_data)
%
% Given a dense d x dN block matrix R whose (dxd)-block elements are
% rotations comprising an optimal solution to the reduced rotation-only
% maximum-likelihood estimation, this function returns a dense (d x n)
% matrix t whose columns are the oorresponding optimal translations

% We have that vec(t) = -(L^pinv V \otimes I_3) vec(R)  (1)
%
% Using the fact that vec(AB) = (B' \otimes I) vec(A),
% we see that (1) is equivalent to
%
% t = - R * V' * L^pinv

Pt = -problem_data.V * R';

% To compute t = P Lw^pinv, observe that 
%
%  t' = (Lw^pinv)' * P' = Lw^pinv * P'   (2) 
%
% since Lw is symmetric.  Now:
%
% Lw = A * Omega * A', and 
%
% =>   Lw^pinv = (A^pinv)' * Omega^{-1} * A^pinv, and 
%
%  A^pinv = Ared' * (Ared * Ared')^{-1} * E
%         = Ared' * (Lred * Lred')^{-1} * E
%         = Ared' * (Lred')^{-1} * Lred^{-1} * E
%
% where E = [I_{n-1} - (1/N)*1_{n-1} 1_{n-1}';   -(1/n)1_{n-1}];

%Compute A1 = E*Pt
Ptop = Pt(1 : problem_data.n - 1, :);
Pbottom = Pt(problem_data.n, :);

A1 = Ptop - (1/problem_data.n)*repmat(sum(Ptop) + Pbottom, problem_data.n-1, 1);

Ptop = Pt(1:problem_data.n - 1, :);
Pbottom = Pt(problem_data.n, :);

%A2 = Lred^{-1} * A1
A2 = problem_data.L \ A1;

%A3 = (Lred')^{-1} * A2
A3 = problem_data.L' \ A2;

%A4 = Ared' * A3
A4 = problem_data.Ared' * A3;

% Middle product in the series
M = problem_data.Omega \ A4;

B1 = problem_data.Ared * M;

%B2 = Lred^{-1} B1
B2 = problem_data.L \ B1;

%B3 = (Lred')^{-1}
B3 = problem_data.L' \ B2;

B4top = B3 - (1/problem_data.n) * repmat(sum(B3), problem_data.n-1, 1);
B4bottom = -(1/problem_data.n)*sum(B3);

B4 = [B4top; B4bottom];

t = B4';

end

