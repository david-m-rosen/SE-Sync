function [B3, B2, B1, B] = construct_B_matrices(measurements)
% function [B3, B2, B1, B] = construct_B_matrices( measurements )
%
% Given a measurement struct of the same format that the SE-Sync method
% expects, this function computes and returns the (sparse) system matrices
% B1, B2, B3, and B defined in the first section of the appendix
% "Reformulating the estimation problem" in the SE-Sync paper

d = length(measurements.t{1});  % Dimension of observations
n = max(max(measurements.edges));  % Number of poses
m = size(measurements.edges, 1);  % Number of measurements

% B3 matrix:
B3_nnz = (d^3 + d^2)*m;
B3_rows = zeros(1, B3_nnz);
B3_cols = zeros(1, B3_nnz);
B3_vals = zeros(1, B3_nnz);

for e = 1:m
    tail = measurements.edges(e,1);
    head = measurements.edges(e,2);
    
    sqkappa = sqrt(measurements.kappa{e});
    Rt = measurements.R{e}';
    
    % Block entries corresponding to the tail of this edge
    
    for r = 1:d
        for c = 1:d
            
            % Block representation of the -sqrt(kappa) *Rt(i,j) * I_d that
            % appears in the Kronecker product
            
            idxs = [(d^3 + d^2)*(e-1) + d^2*(r-1) + d*(c - 1) + 1 : (d^3 + d^2)*(e-1) + d^2*(r-1) + d*(c - 1) + d];
            
            B3_rows(idxs) = [d^2 * (e-1) + d*(r-1) + 1 : d^2 * (e-1) + d*(r-1) + d];
            B3_cols(idxs) = [d^2 * (tail - 1) + d*(c-1) + 1 : d^2 * (tail - 1) + d*(c-1) + d];
            B3_vals(idxs) = -sqkappa * Rt(r,c) * ones(1, d);
        end
    end
    
    % Large d x d diagonal block entry corresponding to the head of this
    % edge
    
    idxs = [d^3 * e + d^2*(e-1) + 1 : (d^3 + d^2)*e];
    
    B3_rows(idxs) = [d^2 *(e-1) + 1 : d^2 * e];
    B3_cols(idxs) = [d^2 *(head -1) + 1 : d^2 * head];
    B3_vals(idxs) = sqkappa * ones(1, d^2);
end

B3 = sparse(B3_rows, B3_cols, B3_vals, d^2*m, d^2*n);

if nargout > 1
    
    % B2 matrix:
    B2_nnz = m*d^2;
    B2_rows = zeros(1, B2_nnz);
    B2_cols = zeros(1, B2_nnz);
    B2_vals = zeros(1, B2_nnz);
    
    for e = 1:m
        tail = measurements.edges(e, 1);
        
        sqtau = sqrt(measurements.tau{e});
        tij = measurements.t{e};
        
        
        for idx = 1:d
            % Block representation of -sqrt(tau)*t_ij(i) * I_d in the Kronecker
            % product:
            
            B2_rows(d^2 * (e-1) + d*(idx-1) + 1 : d^2 * (e-1) + d*idx) = [d*(e-1) + 1 : d*e];
            B2_cols(d^2 * (e-1) + d*(idx-1) + 1 : d^2 * (e-1) + d*idx) = [d^2 * (tail-1) + d*(idx - 1) + 1 : d^2 * (tail-1) + d*idx];
            B2_vals(d^2 * (e-1) + d*(idx-1) + 1 : d^2 * (e-1) + d*idx) = -sqtau * tij(idx)*ones(1,d);
        end
    end
    
    B2 = sparse(B2_rows, B2_cols, B2_vals, d*m, d^2*n);
    
    
    % B1 matrix:
    
    B1_nnz = 2*d*m;
    B1_rows = zeros(1, B1_nnz);
    B1_cols = zeros(1, B1_nnz);
    B1_vals = zeros(1, B1_nnz);
    
    for e = 1:m
        tail = measurements.edges(e,1);
        head = measurements.edges(e,2);
        
        % Block entry corresponding to the tail of this edge
        B1_rows(2*d*(e-1) + 1 : 2*d*e - d) = [d*(e-1) + 1 : d*e];
        B1_cols(2*d*(e-1) + 1 : 2*d*e - d) = [d*(tail - 1) + 1 : d*tail];
        B1_vals(2*d*(e-1) + 1 : 2*d*e - d) = -sqrt(measurements.tau{e})*ones(1, d);
        
        % Block entry corresponding to the head of this edge
        B1_rows(2*d*e - d + 1 : 2*d*e) = [d*(e-1) + 1 : d*e];
        B1_cols(2*d*e - d + 1 : 2*d*e) = [d*(head - 1) + 1 : d*head];
        B1_vals(2*d*e - d + 1 : 2*d*e) = sqrt(measurements.tau{e})*ones(1, d);
    end
    
    B1 = sparse(B1_rows, B1_cols, B1_vals, d*m, d*n);
    
end

if nargout >= 4
    B = [B1, B2;
        sparse(size(B3, 1), size(B1, 2)), B3];
end