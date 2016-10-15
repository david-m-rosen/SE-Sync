function measurements = construct_measurements(filename, is2D)
%function measurements = construct_measurements(filename, is2D)
%
% Given a path 'filename' to a .g2o file, this function reads the file, and
% and constructs the correspoding 'measurements' struct required by SE-Sync


if( nargin < 2 || ~is2D)
    % This is a 3D dataset
    
    % Read in .g2o file
    [g2o_measurements, edges] = readG2oDataset3D(filename, false);
    
    % Store dimensions of the problem
    
    % Dimension of the special Euclidean group
    D = 3;
    % Number of states to estimate
    N = g2o_measurements.nrNodes;
    % Number of available measurements
    M = length(g2o_measurements.between);
    
    % Store the edges in the graph
    measurements.edges = edges;
    
    % Arrays containing the measurements and measurement precisions
    measurements.R = cell(1,M);
    measurements.t = cell(1,M);
    measurements.tau = cell(1,M);
    measurements.kappa = cell(1,M);
    
    for k = 1:M
        measurements.R{k} = g2o_measurements.between(k).R;
        measurements.t{k} = g2o_measurements.between(k).t;
        measurements.tau{k} = g2o_measurements.between(k).Info(1,1);
        measurements.kappa{k} = g2o_measurements.between(k).Info(D+1, D+1);
    end
else
    % This is a 2D dataset
    
    % Read in .g2o file
    [g2o_measurements, readPoses] = DaveReadG2oDataset2D(filename);
    
    D = 2;
    M = size(g2o_measurements, 1);
    N = size(readPoses, 1);
    measurements.edges = g2o_measurements(:, 1:2);
    
    % Arrays containing the measurements and measurement precisions
    measurements.R = cell(1,M);
    measurements.t = cell(1,M);
    measurements.tau = cell(1,M);
    measurements.kappa = cell(1,M);
    
    for k = 1:M
        measurements.t{k} = g2o_measurements(k, 3:4)';
        
        theta = g2o_measurements(k, 5);
        measurements.R{k} = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        % Construct translational measurement information matrix
        I11 = g2o_measurements(k, 6);
        I12 = g2o_measurements(k, 7);
        %I13 = g2o_measurements(k,8);
        I22 = g2o_measurements(k, 9);
        
        Itran = [I11 I12;
                 I12 I22];
             
        measurements.tau{k} = min(eig(Itran));
        measurements.kappa{k} = g2o_measurements(k, 11);  % I33
    end
        
end
    