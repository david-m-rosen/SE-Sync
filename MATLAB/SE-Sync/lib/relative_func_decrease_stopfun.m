function stop = relative_func_decrease_stopfun(manopt_problem, x, info, last, rel_func_decrease_tol)
%function stop = relative_func_decrease_stopfun(manopt_problem, x, info, last, rel_func_decrease_tol)
%
% This function provides an additional stopping criterion for the Manopt
% library: stop whenever the relative decrease in the function value
% between two successive accepted update steps is less than
% rel_func_decrease_tol

% Copyright (C) 2016 by David M. Rosen

this_iterate_accepted = info(last).accepted;

if (~this_iterate_accepted || last == 1)
    stop = false;
else
    % This iterate was an accepted update step, so get the index of the
    % most recent previously-accepted update
    
    accepted_array = [info.accepted];
    previous_accepted_iterate_idx = find(accepted_array(1:last-1), 1, 'last');
    
    if ~isempty(previous_accepted_iterate_idx)
        % Get the function value at the previous accepted iterate
        
        previous_val = info(previous_accepted_iterate_idx).cost;
        current_val = info(last).cost;
        
        if (previous_val - current_val) / previous_val < rel_func_decrease_tol
            stop = true;
        else
            stop = false;
        end
        
    else
        % There have been no previously-accepted update steps
        stop = false;
    end
    
end
end

