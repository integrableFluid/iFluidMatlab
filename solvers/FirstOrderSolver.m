classdef FirstOrderSolver < iFluidSolver
    % Solves GHD Euler-equation using a first order step 

methods (Access = public)
    
    % Superclass constructor
    function obj = FirstOrderSolver(coreObj, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
    end
    
    
end % end public methods


methods (Access = protected)

    % Implementation of abstract functions
    
    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        % =================================================================
        % Purpose : Calculates and stores all quantities needed for the
        %           step() method.
        %           In this case of first order step, nothing is required.
        % Input :   theta_init -- Initial filling function.
        %           u_init     -- Initial position characteristic.
        %           w_init     -- Initial rapidity characteristic.
        %           t_array    -- Array of time steps.
        % Output:   theta      -- Input theta for first step().
        %           u          -- Input u for first step().
        %           w          -- Input w for first step().
        % =================================================================
        theta   = theta_init;
        u       = u_init;
        w       = w_init;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)
        % =================================================================
        % Purpose : Performs a single, first-order Euler step propagating
        %           the filling function theta(t) --> theta(t+dt).
        %           Note, the first order step is already implemented in 
        %           the superclass.
        % Input :   theta_prev -- Filling function at time t.
        %           u_prev     -- Position characteristic at time t.
        %           w_prev     -- Rapidity characteristic at time t.
        %           t          -- Starting time.
        %           dt         -- Length of time step.
        % Output:   theta_next -- Filling function at time t+dt.
        %           u_next     -- Position characteristic at time t+dt.
        %           w_next     -- Rapidity characteristic at time t+dt.
        % =================================================================
        [theta_next, u_next, w_next] = obj.performFirstOrderStep(theta_prev, u_prev, w_prev, t, dt);     
    end
    
end % end protected methods

end % end classdef