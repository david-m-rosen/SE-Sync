function [R, singular_values, determinants] = round_solution(Yopt, problem_data)
%function [R, singular_values, determinants] = round_solution(Yopt, problem_data)
%
% Given an element Yopt in St(d, r)^n, this function rounds Yopt to an
% element of SO(d)^n

% Copyright (C) 2016 by David M. Rosen


r = size(Yopt, 1);

[U, Xi, V] = svd(Yopt, 'econ');
singular_values = diag(Xi)';

Xi_d = Xi(1:problem_data.d, 1:problem_data.d);  %Xi_d is the upper-left dxd submatrix of Xi
V_d = V(:, 1:problem_data.d);  %V_d contains the first d columns of V

R = Xi_d*V_d';


determinants = zeros(1, problem_data.n);

for k = 1:problem_data.n
    determinants(k) = det(R(:, problem_data.d*(k-1) + 1 : problem_data.d*(k-1) + problem_data.d));
end
ng0 = sum(determinants > 0);

reflector = diag([ones(1, problem_data.d - 1), -1]);  % Orthogonal matrix that we can use for reversing the orientations of the orthogonal matrix subblocks of R

if ng0 == 0
    % This solution converged to a reflection of the correct solution
    R = reflector*R;
    determinants = -determinants;
elseif ng0 < problem_data.n
    disp('WARNING: SOLUTION HAS INCONSISTENT ORIENTATIONS!');
    
    % If more than half of the determinants have negative sign, reverse
    % them
    if ng0 < problem_data.n / 2
        determinants = -determinants;
        R = reflector * R;
    end
end

% Finally, project each element of R to SO(d)
for i = 1:problem_data.n
    R(:, problem_data.d * (i-1) + 1 : problem_data.d *i) = project_to_SOd( R(:, problem_data.d * (i-1) + 1 : problem_data.d *i) );
end

end

