classdef AdvectionSolver_BSL_AM4 < AdvectionSolver_BSL
    
    % Backwards Semi-Lagrangian (BSL) solver for the GHD advection equation
    % using an Adams-Moulton scheme of order 4 (AM4). 
    
properties (Access = protected)
    
    % velocity fields of previous steps
    v_m     = [];
    a_m     = [];
    v_mm    = [];
    a_mm    = [];
    v_mmm   = [];
    a_mmm   = [];
    
end % end protected properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSolver_BSL_AM4(model, varargin)        
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
        

        if iscell(theta_init) && length(theta_init) > 3
            % Assume theta_init{end} = theta(t = 0),
            % while theta_init{end-1} = theta(t = -dt), ...
            dt = t_array(2) - t_array(1);
            [obj.v_m, obj.a_m]      = obj.model.calcEffectiveVelocities(theta_init{end-1}, -dt);
            [obj.v_mm, obj.a_mm]    = obj.model.calcEffectiveVelocities(theta_init{end-2}, -2*dt);
            [obj.v_mmm, obj.a_mmm]  = obj.model.calcEffectiveVelocities(theta_init{end-3}, -3*dt);
            
            theta   = theta_init{end};
        else
            [obj.v_m, obj.a_m]      = obj.model.calcEffectiveVelocities(theta_init, 0);
            [obj.v_mm, obj.a_mm]    = obj.model.calcEffectiveVelocities(theta_init, 0);
            [obj.v_mmm, obj.a_mmm]  = obj.model.calcEffectiveVelocities(theta_init, 0);
            
            theta   = theta_init;
        end
    end
      

    function [x_d, r_d, v_n, a_n] = calculateDeparturePoints(obj, theta, t, dt)
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
            v_p         = 4*v_n - 6*obj.v_m + 4*obj.v_mm - obj.v_mmm;
            a_p         = 4*a_n - 6*obj.a_m + 4*obj.a_mm - obj.a_mmm;
        end
        
        x_n         = obj.x_grid - dt*v_n;      % initial point
        r_n         = obj.rapid_grid - dt*a_n;  % initial point
        
        % solve fixed-point problem using Picard iteration    
        error       = 1;
        iter        = 1;
        while error > obj.settings.tol && iter <= obj.settings.max_iter
            v_n_star    = obj.interpPhaseSpace(v_n, r_n, x_n, true);
            a_n_star    = obj.interpPhaseSpace(a_n, r_n, x_n, true);
            
            x_m         = obj.x_grid + 4*(obj.x_grid - x_n) - 2*dt*(2*v_n_star + v_p);
            r_m         = obj.rapid_grid + 4*(obj.rapid_grid - r_n) - 2*dt*(2*a_n_star + a_p);
            
            x_mm        = obj.x_grid + 27*(obj.x_grid - x_n) - 6*dt*(3*v_n_star + 2*v_p);
            r_mm        = obj.rapid_grid + 27*(obj.rapid_grid - r_n) - 6*dt*(3*a_n_star + 2*a_p);
            
            v_m_star    = obj.interpPhaseSpace(obj.v_m, r_m, x_m, true);
            a_m_star    = obj.interpPhaseSpace(obj.a_m, r_m, x_m, true);
            
            v_mm_star   = obj.interpPhaseSpace(obj.v_mm, r_mm, x_mm, true);
            a_mm_star   = obj.interpPhaseSpace(obj.a_mm, r_mm, x_mm, true);
            
            Gx          = obj.x_grid - x_n - dt/24*(9*v_p + 19*v_n_star - 5*v_m_star + v_mm_star);
            Gr          = obj.rapid_grid - r_n - dt/24*(9*a_p + 19*a_n_star - 5*a_m_star + a_mm_star);
            
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
        obj.v_mmm   = obj.v_mm;
        obj.a_mmm   = obj.a_mm;
        obj.v_mm    = obj.v_m;
        obj.a_mm    = obj.a_m;
        obj.v_m     = v;
        obj.a_m     = a;
    end
    
end % end protected methods

end % end classdef