classdef AdvectionSource < handle
    % Superclass for source term of the non-linear GHD advection equation.
    %
    %
properties (Access = protected)
    % Grid lengths
    M                   = []; % number of spatial grid-points
    N                   = []; % number of rapidity grid-points
    Ntypes              = []; % number of quasi-particle types
    
    % Grids (all vectors)
    rapid_grid          = [];
    x_grid              = [];
    type_grid           = [];
    rapid_w             = []; % weights for Gaussian quadrature

    model               = []; % iFluidCore object specifying the model
    SolverObj           = []; % AvectionSolver_BSL object for interpolation
    settings            = []; % struct specifying settings
    
    S_prev              = []; % main source term of previous step
    A_prev              = []; % auxiliary source term of previous step
    
end % end protected properties



methods (Abstract, Access = public)
    
    % Abstract methods must be implemented in extending class!
    
    [fill, u, w]   = initialize(obj, fill_init, u_init, w_init, t_array)
    [S, A]          = calcSourceTerm(obj, fill, u, w, t);
    [u_next, w_next]= auxStep(obj, u, w, A, t, dt); 

end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSource(model, varargin)        
        % iFluidSolver requires an iFluidCore object for calculating
        % velocities of quasiparticles etc.
        assert( isa( model, 'iFluidCore' ) )
        
        % Copy grids from iFluidCore object
        [x_grid, rapid_grid, type_grid, rapid_w] = model.getGrids();

        obj.x_grid      = x_grid;
        obj.rapid_grid  = rapid_grid;
        obj.type_grid   = type_grid;
        obj.rapid_w     = rapid_w;
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = length(type_grid);
        
        obj.model       = model;     
        
        
        % Parse inputs from varargin
        parser = inputParser;
        parser.KeepUnmatched = true;

        macthStrings        = @(x,str) any(validatestring(x,str));
        default_scheme      = 'SETTLS';
        expected_scheme     = {'Crank-Nicolson','SETTLS','midpoint','endpoint'};
        addParameter(parser,'propagation_scheme',default_scheme, @(x) macthStrings(x,expected_scheme))
        
        addParameter(parser, 'tol', 1e-4, @(x) x > 0);
        addParameter(parser, 'max_iter', 50, @(x) isscalar(x) & floor(x) == x & x >=0);
        
        parse(parser,varargin{:});
        
        % set settings equal to parser results
        obj.settings  = parser.Results;
    end
    

    function attachSolver(obj, solver)
        obj.SolverObj = solver;
    end
    
    
    function [fill_next, u_next, w_next] = step(obj, fill, x_d, r_d, u, w, t, dt)
        % =================================================================
        % Purpose : Propagates filling function fill one timestep dt 
        %           according to the advection source term.
        % Input :   fill   -- Filling function at time t.
        %           x_d     -- Departure points for time-step dt.
        %           r_d     -- Departure points for time-step dt.
        %           u       -- Auxiliary variable.
        %           w       -- Auxiliary variable.
        %           t       -- Time
        %           dt      -- Time-step
        % Output:   fill_next -- Propagated filling.
        %           u_next     -- Propagated u.
        %           w_next     -- Propagated w.
        % =================================================================
        
        
        switch obj.settings.propagation_scheme
            
        case 'Crank-Nicolson'
            
            % interpolate filling function to departure points
            fill_s     = obj.SolverObj.interpPhaseSpace(fill, r_d, x_d);
            
            % calculate source term at start-point
            [S, A_s]    = obj.calcSourceTerm(fill, u, w, t);
            S_s         = obj.SolverObj.interpPhaseSpace(S, r_d, x_d);
%             [S_d, A_d]  = obj.calcSourceTerm(fill_d, u, w, t);

            % solve fixed-point problem (G -> 0) using Picard iteration
            fill_e      = fill_s;  % initial point
            u_a         = u;
            w_a         = w;
            err       	= 1;
            iter        = 1;
            while err > obj.settings.tol && iter <= obj.settings.max_iter
                % calculate source term at end-point
                [S_e, A_e] = obj.calcSourceTerm(fill_e, u_a, w_a, t+dt);

                % setup fixed-point problem 
                G           = fill_s - fill_e + dt/2*( S_s + S_e );
                fill_e     = fill_e + G;
                
                % update auxiliary functions using the mid-point of the
                % auxiliary source term A
                [u_a, w_a]  = obj.auxStep(u, w, 0.5*(A_s + A_e), t, dt);
                
                % calculate error for termination condition
                err         = double(sum(G.^2,'all'));
                iter        = iter + 1;
            end

            % return filling and auxiliary functions at arrival point
            fill_next  = fill_e;
            u_next      = u_a;
            w_next      = w_a;  

        case 'SETTLS'
            
            % interpolate filling function to departure points
            fill_s     = obj.SolverObj.interpPhaseSpace(fill, r_d, x_d);
            
            % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill, u, w, t);
            
            % approximate source term at time t+dt and interpolate to
            % departure points
            S_end_d     = obj.SolverObj.interpPhaseSpace((2*S-obj.S_prev), r_d, x_d);
            
            % approximate auxiliary source term A at time t+dt/2
            A_mid       = 1.5*A - 0.5*obj.A_prev;

            % update filling and auxiliary functions
            fill_next      = fill_s + dt/2*( S + S_end_d );
            [u_next, w_next]= obj.auxStep(u, w, A_mid, t, dt);        

            % store results for next iteration
            obj.S_prev      = S;
            obj.A_prev      = A;
            
        case 'midpoint'
            
            % interpolate filling function to departure points
            fill_s     = obj.SolverObj.interpPhaseSpace(fill, r_d, x_d);
            
             % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill, u, w, t);
            
            % approximate source terms A and A at time t+dt/2
            S_mid           = 1.5*S - 0.5*obj.S_prev;
            A_mid           = 1.5*A - 0.5*obj.A_prev;
            S_mid_d         = obj.SolverObj.interpPhaseSpace(S_mid, r_d, x_d);

            % update filling and auxiliary functions
            fill_next      = fill_s + dt/2*( S_mid + S_mid_d );
            [u_next, w_next]= obj.auxStep(u, w, A_mid, t, dt);        

            % store results for next iteration
            obj.S_prev      = S;
            obj.A_prev      = A;            
            
        case 'endpoint'
            
            % interpolate filling function to departure points
            fill_s     = obj.SolverObj.interpPhaseSpace(fill, r_d, x_d);
            
            % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill_s, u, w, t);
        
            % update filling and auxiliary functions
            fill_next      = fill_s + dt*S;
            [u_next, w_next]= obj.auxStep(u, w, A, t, dt); 

        otherwise
            error('Detected source term, but no valid hybrid scheme is selected.')
        end
        
    end



end % end public methods

end % end classdef