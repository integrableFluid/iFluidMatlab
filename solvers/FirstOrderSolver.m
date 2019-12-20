classdef FirstOrderSolver < iFluidSolver


methods (Access = public)
    
    % Superclass constructor
    function obj = FirstOrderSolver(coreObj, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
    end
    
    
end % end public methods


methods (Access = protected)

    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        % No special initialization required for first order step
        theta   = theta_init;
        u       = u_init;
        w       = w_init;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)
        % Step function is simply the already implmented first order step
        [theta_next, u_next, w_next] = obj.performFirstOrderStep(theta_prev, u_prev, w_prev, t, dt);     
    end
    
end % end protected methods

end % end classdef