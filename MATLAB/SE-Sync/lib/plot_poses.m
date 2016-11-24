function f = plot_poses(t_hat, Rhat, edges, lc_linestyle, lc_alpha)
%function f = plot_poses(t_hat, Rhat, edges, lc_linestyle, lc_alpha)
%
% Given translational and rotation state estimates returned by SE_Sync,
% this function plots the corresponding solution, applying the gauge
% symmetry that maps the first pose to the identity element of SE(d).  The
% 'edges' argument is optional; if supplied, it will plot loop closures;
% otherwise only odometric links are shown.  lc_linestyle is an optional
% string specifying the line style to be used for plotting the loop closure
% edges (default: '-b'); similarly, lc_alpha is an optional alpha value for
% the loop closure edges (default: alpha = 1.0)

if(~exist('lc_linestyle'))
    lc_linestyle = '-b';
end
if(~exist('lc_alpha'))
    lc_alpha = 1.0;
end

D = size(t_hat, 1);  % Get dimension of solutions, should be either 2 or 3

% 'Unrotate' these vectors by premultiplying by the inverse of the first
% orientation estimate
t_hat_rotated = Rhat(1:D, 1:D)' * t_hat;
% Translate the resulting vectors to the origin
t_hat_anchored = t_hat_rotated - repmat(t_hat_rotated(:, 1), 1, size(t_hat_rotated, 2));

x = t_hat_anchored(1,:);
y = t_hat_anchored(2, :);

if D == 3
    z = t_hat_anchored(3, :);
    
    % Plot odometric links
    f = figure();
    plot3(x, y, z, '-b');
    axis equal;
    hold on;
    
    if ( nargin > 2)
        
        for k = 1:size(edges, 1)
            id1 = edges(k, 1);
            id2 = edges(k, 2);
            
            if abs(id1 - id2) > 1
                % This is a loop closure measurement
                lc_plot = plot3(t_hat_anchored(1, [id1 id2]), t_hat_anchored(2, [id1 id2]), t_hat_anchored(3, [id1 id2]), lc_linestyle, 'Linewidth', 1);
                lc_plot.Color(4) = lc_alpha;  % Set transparency of loop closure edges
            end
        end
    end
elseif D == 2
    
    % Plot odometric links
    f = figure();
    plot(x, y, '-b');
    axis equal;
    hold on;
    
    if ( nargin > 2)
        
        for k = 1:size(edges, 1)
            id1 = edges(k, 1);
            id2 = edges(k, 2);
            
            if abs(id1 - id2) > 1
                % This is a loop closure measurement
                lc_plot = plot(t_hat_anchored(1, [id1 id2]), t_hat_anchored(2, [id1 id2]), lc_linestyle);
                lc_plot.Color(4) = lc_alpha;  % Set transparency of loop closure edges
            end
        end
    end
    
end
end

