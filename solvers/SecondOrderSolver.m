classdef SecondOrderSolver < iFluidSolver

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
        dt      = t_array(2) - t_array(1);
        ddt     = dt/2/10;
        theta   = theta_init;
        u       = u_init;
        w       = w_init;

        % Calculate first theta_mid at t = dt/10/2 using first order
        % step, then use that to calculate the actual theta_mid at
        % t = dt/2 using second order steps
        
        obj.theta_mid = obj.performFirstOrderStep(theta_init, u_init, w_init, 0, ddt);
        theta_temp = theta;

        for i = 1:10
            t           = (i-1)*ddt;
            theta_temp  = obj.step(theta_temp, u_init, w_init, t, ddt);
        end

        obj.theta_mid = theta_temp;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_next, t, dt)
        % Calculate theta^[n+1] using the filling (theta) at the midpoint
        % between theta^[n] and theta^[n+1].
        [theta_next, u_next, w_next] = step2(obj, obj.theta_mid, theta_prev, u_prev, w_next, t, dt);
        
        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % Doesnt output u_mid, as it is not needed.
        obj.theta_mid  = step2(obj, theta_next, obj.theta_mid, u_prev, w_next, t+dt/2, dt);
        
        
        function [theta_next, u_next, w_next] = step2(obj, theta_mid, theta_prev, u_prev, w_prev, t, dt)
            % Estimate x' and rapid' using midpoint filling
            [v_eff, a_eff] = obj.coreObj.calcEffectiveVelocities(theta_mid, t+dt/2, obj.x_grid, obj.rapid_grid, obj.type_grid);

            x_mid   = obj.x_grid - 0.5*dt*v_eff; % (1xNxM)
            r_mid   = obj.rapid_grid - 0.5*dt*a_eff; % (1xNxM)

            % Interpolate v_eff and a_eff to midpoint coordinates x' and rapid' 
            v_mid   = obj.interpPhaseSpace( v_eff, r_mid, x_mid, true ); % always extrapolate v_eff
            a_mid   = obj.interpPhaseSpace( a_eff, r_mid, x_mid, true );
            
            % From midpoint velocity calculate fullstep x and rapid translation
            x_back  = obj.x_grid - dt*v_mid;
            r_back  = obj.rapid_grid - dt*a_mid;

            % Use interpolation to find theta_prev at x_back, r_back and
            % assign values to theta_next.
            theta_next = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
            u_next     = obj.interpPhaseSpace(u_prev, r_back, x_back, true); % always extrapolate u
            w_next     = obj.interpPhaseSpace(w_prev, r_back, x_back, true); % always extrapolate u
        end % end nested function
    end
    
end % end protected methods

end % end classdef