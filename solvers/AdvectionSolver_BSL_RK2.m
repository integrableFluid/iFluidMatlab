classdef AdvectionSolver_BSL_RK2 < AdvectionSolver_BSL
    
    % Backwards Semi-Lagrangian (BSL) solver for the GHD advection equation
    % using a Runge-Kutte scheme of order 2 (RK2).    
        

methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSolver_BSL_RK2(model, varargin)        
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

        % calculate derivatives of velocities and use to estimate
        % velocities at temporal midpoint
        V           = obj.calcVelocityDerivatives(obj.settings.deriv_order, theta, t);

        veff_mid    = V{1,1};
        aeff_mid    = V{2,1};
        for i = 1:obj.settings.deriv_order
            veff_mid    = veff_mid + (dt/2)^i * V{1,i+1}/factorial(i);
            aeff_mid    = aeff_mid + (dt/2)^i * V{2,i+1}/factorial(i);
        end
        
        % not used
        v_n     = 0;
        a_n     = 0;
        
        
        if obj.settings.implicit
            % solve fixed-point problem for characteristics using Picard
            % iteration.      
            x_d         = obj.x_grid - dt*V{1,1};     % initial point
            r_d         = obj.rapid_grid - dt*V{2,1}; % initial point
            error       = 1;
            iter        = 1;
            while error > obj.settings.tol && iter <= obj.settings.max_iter
                veff_d      = obj.interpPhaseSpace(veff_mid, (obj.rapid_grid+r_d)/2, (obj.x_grid+x_d)/2, true, false);
                aeff_d      = obj.interpPhaseSpace(aeff_mid, (obj.rapid_grid+r_d)/2, (obj.x_grid+x_d)/2, true, false);

                Gx          = obj.x_grid - x_d - dt*veff_d;
                Gr          = obj.rapid_grid - r_d - dt*aeff_d;

                x_d         = x_d + Gx;
                r_d         = r_d + Gr;    

                error       = double(sum(Gx.^2,'all') + sum(Gr.^2,'all'));
                iter        = iter + 1;
            end
        else
            % solve using explicit scheme
            x_star      = obj.x_grid - dt/2*V{1,1};
            r_star      = obj.rapid_grid - dt/2*V{2,1};

            veff_star   = obj.interpPhaseSpace(veff_mid, r_star, x_star, true, false);
            aeff_star   = obj.interpPhaseSpace(aeff_mid, r_star, x_star, true, false);
        
            x_d         = obj.x_grid - dt*veff_star;
            r_d         = obj.rapid_grid - dt*aeff_star;
        end
        
    end
    
    
    function storeVelocityFields(obj, theta, x_d, r_d, v, a)
        % RK solvers do not rely on velocity fields of previous steps
    end
    
end % end protected methods

end % end classdef