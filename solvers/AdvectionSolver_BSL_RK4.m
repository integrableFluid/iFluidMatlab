classdef AdvectionSolver_BSL_RK4 < AdvectionSolver_BSL
    
    % Backwards Semi-Lagrangian (BSL) solver for the GHD advection equation
    % using a Runge-Kutte scheme of order 4 (RK$).    
        

methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSolver_BSL_RK4(model, varargin)        
        obj = obj@AdvectionSolver_BSL(model, varargin{:});
        
        % Parse inputs from varargin
        parser = inputParser;
        parser.KeepUnmatched = true;

        addParameter(parser, 'implicit', false, @(x) islogical(x));

        parse(parser,varargin{:});
        
        % add parser results to settings
        names = [fieldnames(obj.settings); fieldnames(parser.Results)];
        obj.settings = cell2struct([struct2cell(obj.settings); struct2cell(parser.Results)], names, 1);
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
      

    function [x_d, r_d, v_n, a_n] = calculateDeparturePoints(obj, theta, t, dt)
        % Note, for interpolation of velocities, always extrapolate (true)
        % but ignore boundary conditions (false).


        % not used
        v_n     = 0;
        a_n     = 0;        
        
        % calculate derivatives of velocities and use to estimate
        % velocities at temporal midpoint
        V           = obj.calcVelocityDerivatives(obj.settings.deriv_order, theta, t);

        veff_mid    = V{1,1};
        aeff_mid    = V{2,1};
        veff_end    = V{1,1};
        aeff_end    = V{2,1};
        for i = 1:obj.settings.deriv_order
            veff_mid    = veff_mid + (dt/2)^i * V{1,i+1}/factorial(i);
            aeff_mid    = aeff_mid + (dt/2)^i * V{2,i+1}/factorial(i);
            
            veff_end    = veff_end + dt^i * V{1,i+1}/factorial(i);
            aeff_end    = aeff_end + dt^i * V{2,i+1}/factorial(i);
        end
        
        
        if obj.settings.implicit
            % solve fixed-point problem for characteristics using Picard
            % iteration.                 
            x_d         = obj.x_grid - dt*(V{1,1} + veff_end)/2;     % initial point
            r_d         = obj.rapid_grid - dt*(V{2,1} + aeff_end)/2; % initial point
            error       = 1;
            iter        = 1;
            while error > obj.settings.tol && iter <= obj.settings.max_iter
                % departure point velocities
                v1          = obj.interpPhaseSpace(V{1,1}, r_d, x_d, true, false);
                a1          = obj.interpPhaseSpace(V{2,1}, r_d, x_d, true, false);
                
                % first midpoint velocities
                v2          = obj.interpPhaseSpace(veff_mid, r_d+dt/2*a1, x_d+dt/2*v1, true, false);
                a2          = obj.interpPhaseSpace(aeff_mid, r_d+dt/2*a1, x_d+dt/2*v1, true, false);
                
                % second midpoint velocities
                v3          = obj.interpPhaseSpace(veff_mid, r_d+dt/2*a2, x_d+dt/2*v2, true, false);
                a3          = obj.interpPhaseSpace(aeff_mid, r_d+dt/2*a2, x_d+dt/2*v2, true, false);

                % arrival point velocities
                v4          = veff_end;
                a4          = aeff_end; 
                
                Gx          = obj.x_grid - x_d - dt/6*(v1 + 2*v2 + 2*v3 + v4);
                Gr          = obj.rapid_grid - r_d - dt/6*(a1 + 2*a2 + 2*a3 + a4);

                x_d         = x_d + Gx;
                r_d         = r_d + Gr;    

                error       = double(sum(Gx.^2,'all') + sum(Gr.^2,'all'));
                iter        = iter + 1;
            end

        else
            % solve using explicit scheme
            
            % first midpoint estimate
            x_star1     = obj.x_grid - dt/2*V{1,1};
            r_star1     = obj.rapid_grid - dt/2*V{2,1};
            veff_star1  = obj.interpPhaseSpace(veff_mid, r_star1, x_star1, true, false);
            aeff_star1  = obj.interpPhaseSpace(aeff_mid, r_star1, x_star1, true, false);
        
            % second midpoint estimate
            x_star2     = obj.x_grid - dt/2*veff_star1;
            r_star2     = obj.rapid_grid - dt/2*aeff_star1;
            veff_star2  = obj.interpPhaseSpace(veff_mid, r_star2, x_star2, true, false);
            aeff_star2  = obj.interpPhaseSpace(aeff_mid, r_star2, x_star2, true, false);
            
            % endpoint estimate
            x_star3     = obj.x_grid - dt*veff_star2;
            r_star3     = obj.rapid_grid - dt*aeff_star2;
            veff_star3  = obj.interpPhaseSpace(veff_end, r_star3, x_star3, true, false);
            aeff_star3  = obj.interpPhaseSpace(aeff_end, r_star3, x_star3, true, false);
            
            % calculate departure points
            x_d         = obj.x_grid - dt/6*( V{1,1} + 2*veff_star1 + 2*veff_star2 + veff_star3 );
            r_d         = obj.rapid_grid - dt/6*( V{2,1} + 2*aeff_star1 + 2*aeff_star2 + aeff_star3 );
        end
        
    end
    
    
    function storeVelocityFields(obj, theta, x_d, r_d, v, a)
        % RK solvers do not rely on velocity fields of previous steps
    end
    
end % end protected methods

end % end classdef