classdef AdvectionSolver_BSL_AM2 < AdvectionSolver_BSL
    
    % Backwards Semi-Lagrangian (BSL) solver for the GHD advection equation
    % using an Adams-Moulton scheme of order 2 (AM2). 
    
properties (Access = protected)
    
    % velocity fields of previous steps
    v_m     = [];
    a_m     = [];
    
end % end protected properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSolver_BSL_AM2(model, varargin)        
        obj = obj@AdvectionSolver_BSL(model, varargin{:});
        
        % Parse inputs from varargin
        parser = inputParser;
        parser.KeepUnmatched = true;

        addParameter(parser, 'extrap_velocity', true, @(x) islogical(x));

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
        u       = u_init;
        w       = w_init;
        

        if iscell(theta_init) && length(theta_init) > 1
            % Assume theta_init{end} = theta(t = 0),
            % while theta_init{end-1} = theta(t = -dt), ...
            dt = t_array(2) - t_array(1);
            [obj.v_m, obj.a_m]      = obj.model.calcEffectiveVelocities(theta_init{end-1}, -dt);
            
            theta   = theta_init{end};
        else
            [obj.v_m, obj.a_m]      = obj.model.calcEffectiveVelocities(theta_init, 0);
            
            theta   = theta_init;
        end
    end
      

    function [x_d, r_d, v_n, a_n] = calculateDeparturePoints(obj, theta, t, dt)
        % Note, for interpolation of velocities, always extrapolate (true)
        % but ignore boundary conditions (false).

        if ~obj.settings.extrap_velocity
            % use analytic derivative to estimate velocities at next step
            V      = obj.calcVelocityDerivatives(obj.settings.deriv_order, theta, t);

            v_n    = V{1,1};
            a_n    = V{2,1};
            
            v_p    = v_n;
            a_p    = a_n;
            for i = 1:obj.settings.deriv_order
                v_p    = v_p + (dt)^i * V{1,i+1}/factorial(i);
                a_p    = a_p + (dt)^i * V{2,i+1}/factorial(i);
            end
        else
            % use extrapolation to estimate velocities at next step
            [v_n, a_n]  = obj.model.calcEffectiveVelocities(theta, t);
            v_p         = 2*v_n - obj.v_m;
            a_p         = 2*a_n - obj.a_m;
        end
        
        x_n         = obj.x_grid - dt*v_n;      % initial point
        r_n         = obj.rapid_grid - dt*a_n;  % initial point
        
        % solve fixed-point problem using Picard iteration    
        error       = 1;
        iter        = 1;
        while error > obj.settings.tol && iter <= obj.settings.max_iter
            v_n_star    = obj.interpPhaseSpace(v_n, r_n, x_n, true, false);
            a_n_star    = obj.interpPhaseSpace(a_n, r_n, x_n, true, false);
            
            Gx          = obj.x_grid - x_n - dt/2*(v_p + v_n_star);
            Gr          = obj.rapid_grid - r_n - dt/2*(a_p + a_n_star);
            
            x_n          = x_n + Gx;
            r_n          = r_n + Gr;    
            
            error       = double(sum(Gx.^2,'all') + sum(Gr.^2,'all'));
            iter        = iter + 1;
        end
        
        x_d         = x_n;
        r_d         = r_n;
        
    end
    
    
    function storeVelocityFields(obj, theta, x_d, r_d, v, a)
        
        % store previous velocities
        obj.v_m     = v;
        obj.a_m     = a;
    end
    
end % end protected methods

end % end classdef