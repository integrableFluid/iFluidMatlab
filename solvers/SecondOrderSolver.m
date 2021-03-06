classdef SecondOrderSolver < iFluidSolver
    % Solves GHD Euler-equation using a second order step 
    
properties (Access = protected)
    theta_mid = []; % Midpoint filling required for taking 2nd order step

end % end protected properties


methods (Access = public)
    
    % Superclass constructor
    function obj = SecondOrderSolver(coreObj, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
    end
    
    
end % end public methods


methods (Access = protected)

    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        % =================================================================
        % Purpose : Calculates and stores all quantities needed for the
        %           step() method.
        %           In this the filling at dt/2 is needed to start the
        %            propagation.
        % Input :   theta_init -- Initial filling function.
        %           u_init     -- Initial position characteristic.
        %           w_init     -- Initial rapidity characteristic.
        %           t_array    -- Array of time steps.
        % Output:   theta      -- Input theta for first step().
        %           u          -- Input u for first step().
        %           w          -- Input w for first step().
        % =================================================================
        dt      = t_array(2) - t_array(1);
        ddt     = dt/2/10;
        theta   = theta_init;
        u       = u_init;
        w       = w_init;

        % Calculate first theta_mid at t = dt/10/2 using first order
        % step, then use that to calculate the actual theta_mid at
        % t = dt/2 using second order steps. 
        obj.theta_mid = obj.performFirstOrderStep(theta_init, u_init, w_init, 0, ddt/2);
        theta_temp = theta;

        for i = 1:10
            t           = (i-1)*ddt;
            theta_temp  = obj.step(theta_temp, u_init, w_init, t, ddt);
        end

        obj.theta_mid = theta_temp;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_next, t, dt)
        % =================================================================
        % Purpose : Performs a single, first-order Euler step propagating
        %           the filling function theta(t) --> theta(t+dt).
        % Input :   theta_prev -- Filling function at time t.
        %           u_prev     -- Position characteristic at time t.
        %           w_prev     -- Rapidity characteristic at time t.
        %           t          -- Starting time.
        %           dt         -- Length of time step.
        % Output:   theta_next -- Filling function at time t+dt.
        %           u_next     -- Position characteristic at time t+dt.
        %           w_next     -- Rapidity characteristic at time t+dt.
        % =================================================================
        
        % First, calculate theta(t+dt) using the midpoint filling theta(t+dt/2)
        [theta_next, u_next, w_next] = step2(obj, obj.theta_mid, theta_prev, u_prev, w_next, t, dt);
        
        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % I.e. calculate theta(t+dt+dt/2) using theta(t+dt) as midpoint.
        obj.theta_mid  = step2(obj, theta_next, obj.theta_mid, u_prev, w_next, t+dt/2, dt);
        
        
        function [theta_next, u_next, w_next] = step2(obj, theta_mid, theta_prev, u_prev, w_prev, t, dt)
            % Estimate x' and rapid' using midpoint filling
            [v_eff, a_eff] = obj.coreObj.calcEffectiveVelocities(theta_mid, t+dt/2, obj.x_grid, obj.rapid_grid, obj.type_grid);

            x_mid   = obj.x_grid - 0.5*dt*v_eff;
            r_mid   = obj.rapid_grid - 0.5*dt*a_eff; 

            % Interpolate v_eff and a_eff to midpoint coordinates x' and rapid' 
            v_mid   = obj.interpPhaseSpace( v_eff, r_mid, x_mid, true ); % always extrapolate v_eff
            a_mid   = obj.interpPhaseSpace( a_eff, r_mid, x_mid, true );
            
            % From midpoint velocity calculate fullstep x and rapid translation
            x_back  = obj.x_grid - dt*v_mid;
            r_back  = obj.rapid_grid - dt*a_mid;

            % Use interpolation to find theta_prev at x_back, r_back and
            % assign values to theta_next.
            theta_next = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
            
            if obj.calcCharac
                u_next  = obj.interpPhaseSpace(u_prev, r_back, x_back, true); % always extrapolate u
                w_next  = obj.interpPhaseSpace(w_prev, r_back, x_back, true); % always extrapolate w
            else % calculate trajectories
                v_int = obj.interpPhaseSpace(v_eff, w_prev, u_prev, true);
                a_int = obj.interpPhaseSpace(a_eff, w_prev, u_prev, true);
                u_next  = u_prev + dt*v_int;
                w_next  = w_prev + dt*a_int;         
            end
        end % end nested function
    end
    
end % end protected methods

end % end classdef