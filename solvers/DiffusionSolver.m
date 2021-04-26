classdef DiffusionSolver < iFluidSolver
    
% Solves GHD equation with diffusion term using a split-step propagation
% scheme.
    
properties (Access = protected)
    theta_mid = []; % Midpoint filling required for taking 2nd order step
    
    F_prev = []; 
    F_mid = [];
    
    
end % end protected properties


methods (Access = public)
    
    % Superclass constructor
    function obj = DiffusionSolver(coreObj, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
        
        obj.calcCharac = false;
                
    end
    
    
    function [F, D] = calcDiffusion(obj, rho, rhoS, t)  
        % =================================================================
        % Purpose : Calculates diffusion term of the GHD equation.
        % Input :   rho        -- Quasiparticle density at time t.
        %           rhoS       -- Density of states at time t.
        %           t          -- Time (scalar).
        % Output:   F          -- Diffusion term.
        %           D          -- Diffusion kernel.
        % =================================================================
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        
        theta       = rho./rhoS;
        v_eff       = obj.coreObj.calcEffectiveVelocities(theta, t, obj.x_grid, obj.rapid_grid, obj.type_grid);  
        I           = fluidcell.eye(obj.N, obj.Ntypes)./obj.rapid_w;
        T           = -1/(2*pi)*obj.coreObj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);
        T_dr        = obj.coreObj.applyDressing(T, theta, t);     
        f           = obj.coreObj.getStatFactor( theta );
        
        % Calculate diffusion kernel       
        W           = rho.*f.*T_dr.^2 .* abs(v_eff - v_eff.t());
        w           = sum( transpose(W.*obj.rapid_w) , 4);
        R           = (I - theta.*T)./rhoS;  
        DT          = rhoS.^(-2).*(I.*w - W);
        
        D           = inv(R)*(obj.rapid_w.* ( DT*(R.*obj.rapid_w)));
        dx_rho      = gradient(rho, 'x', obj.x_grid);
        F           = 0.5*gradient( D*(dx_rho), 'x', obj.x_grid );
        
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
        dt          = t_array(2) - t_array(1);
        ddt         = dt/2/10;
        theta       = theta_init;
        u           = u_init;          
        w           = w_init;
        
        
        % Set up collision integrals for propagation
        [rho, rhoS] = obj.coreObj.transform2rho(theta_init, 0);
        F           = obj.calcDiffusion( rho, rhoS, 0);
        obj.F_prev  = F;
        obj.F_mid   = F;

        % Calculate first theta_mid at t = dt/10/2 using first order
        % step, then use that to calculate the actual theta_mid at
        % t = dt/2 using second order steps.
        obj.theta_mid = theta_init + ddt/2*F;
        obj.theta_mid = obj.performFirstOrderStep(obj.theta_mid, u_init, w_init, 0, ddt/2);

        
        theta_temp  = theta;
        
        for i = 1:10
            t          = (i-1)*ddt;
            theta_temp = obj.step(theta_temp, u_init, w_init, t, ddt);
        end

        obj.theta_mid= theta_temp;
        obj.F_mid    = obj.F_prev;
        obj.F_prev   = F;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_next, t, dt)
        % =================================================================
        % Purpose : Performs a split step propagating
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
        
        u_next      = 0;
        w_next      = 0;
        
        % Calculate collision integral
        [rho, rhoS] = obj.coreObj.transform2rho(theta_prev, t);
        F           = obj.calcDiffusion( rho, rhoS, t);
        
        % First part of split step: solve collision equation
        rho         = rho + 0.5*dt*(3*F - obj.F_prev);
        obj.F_prev  = F;
        theta_prev  = obj.coreObj.transform2theta(rho, t);
        
        % Second part of split step: solve propagation equation
        theta_next = step2(obj, obj.theta_mid, theta_prev, t, dt);


        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % I.e. calculate theta(t+dt+dt/2) using theta(t+dt) as midpoint.
        [rho, rhoS] = obj.coreObj.transform2rho(obj.theta_mid, t+dt/2);
        F           = obj.calcDiffusion( rho, rhoS, t+dt/2);
        rho         = rho + 0.5*dt*(3*F - obj.F_mid);
        obj.F_mid   = F;
        
        obj.theta_mid = obj.coreObj.transform2theta(rho, t+dt/2);
        obj.theta_mid = step2(obj, theta_next, obj.theta_mid, t+dt/2, dt);
        

        function theta_next = step2(obj, theta_mid, theta_prev, t, dt)
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

        end % end nested function
    end
    
    
end % end protected methods


methods (Access = private)
   
    
end % end private methods

end % end classdef