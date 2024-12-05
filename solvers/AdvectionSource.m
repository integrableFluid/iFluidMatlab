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
        
        addParameter(parser, 'extrapolate', false, @(x) islogical(x));
        addParameter(parser, 'periodic_rapid', false, @(x) islogical(x));
        addParameter(parser, 'periodic_BC', false, @(x) islogical(x));
        addParameter(parser, 'reflective_BC', false, @(x) islogical(x));
        addParameter(parser, 'tol', 1e-4, @(x) x > 0);
        addParameter(parser, 'max_iter', 50, @(x) isscalar(x) & floor(x) == x & x >=0);
        
        parse(parser,varargin{:});
        
        % set settings equal to parser results
        obj.settings  = parser.Results;
    end
    
    
    function copySettings(obj, settings_in)
        % For each field in struct settings_in, if also in obj.settings,
        % copy it to obj.settings.
        % Used to have matching settings between advection solver and
        % source term.
        
        fn = fieldnames(settings_in);
        for k = 1:numel(fn) % loop over all fields in settings_in
            if isfield(obj.settings, fn{k}) % check if field is also in settings of source
                obj.settings.(fn{k}) = settings_in.(fn{k}); % copy value
            end
        end
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
            fill_s     = obj.interpPhaseSpace(fill, r_d, x_d, obj.settings.extrapolate);
            
            % calculate source term at start-point
            [S, A_s]    = obj.calcSourceTerm(fill, u, w, t);
            S_s         = obj.interpPhaseSpace(S, r_d, x_d, obj.settings.extrapolate);
%             [S_d, A_d]  = obj.calcSourceTerm(fill_d, u, w, t);

            % solve fixed-point problem (G -> 0) using Picard iteration
            fill_e     = fill_s;  % initial point
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
            fill_s     = obj.interpPhaseSpace(fill, r_d, x_d, obj.settings.extrapolate);
            
            % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill, u, w, t);
            
            % approximate source term at time t+dt and interpolate to
            % departure points
            S_end_d     = obj.interpPhaseSpace((2*S-obj.S_prev), r_d, x_d, obj.settings.extrapolate);
            
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
            fill_s     = obj.interpPhaseSpace(fill, r_d, x_d, obj.settings.extrapolate);
            
             % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill, u, w, t);
            
            % approximate source terms A and A at time t+dt/2
            S_mid           = 1.5*S - 0.5*obj.S_prev;
            A_mid           = 1.5*A - 0.5*obj.A_prev;
            S_mid_d         = obj.interpPhaseSpace(S_mid, r_d, x_d, obj.settings.extrapolate);

            % update filling and auxiliary functions
            fill_next      = fill_s + dt/2*( S_mid + S_mid_d );
            [u_next, w_next]= obj.auxStep(u, w, A_mid, t, dt);        

            % store results for next iteration
            obj.S_prev      = S;
            obj.A_prev      = A;            
            
        case 'endpoint'
            
            % interpolate filling function to departure points
            fill_s     = obj.interpPhaseSpace(fill, r_d, x_d, obj.settings.extrapolate);
            
            % calculate source term on grid at time t
            [S, A]      = obj.calcSourceTerm(fill_s, u, w, t);
        
            % update filling and auxiliary functions
            fill_next      = fill_s + dt*S;
            [u_next, w_next]= obj.auxStep(u, w, A, t, dt); 

        otherwise
            error('Detected source term, but no valid hybrid scheme is selected.')
        end
        
    end

    
    function tensor_int = interpPhaseSpace(obj, tensor_grid, rapid_int, x_int, extrapolate)
        % =================================================================
        % Purpose : Interpolates an iFluidTensor defined on the grids 
        %           stored in the object to new coordinates.
        %           This function exists because MATLAB has different
        %           syntax between interp1 and interp2in terms of 
        %           extrapolation...
        % Input :   tensor_grid -- iFluidTensor defined on rapid_grid,
        %                          x_grid, and type_grid.
        %           rapid_int   -- Rapidity values to interpolate to
        %                           (should be iFluidTensor sized
        %                            [N, M, Ntypes] )
        %           x_int       -- Spatial values to interpolate to
        %                           (should be iFluidTensor sized
        %                            [N, M, Ntypes] )
        %           extrapFlag  -- if true, enable extrapolations
        %                          if false, all extrap. values are zero
        % Output:   tensor_int -- iFluidTensor interpolated to input grids.
        % =================================================================
        

        % interpolation method
        method = 'spline';
        if GPU_mode_on()
            method = 'cubic';
        end

        % Cast to matrix form
        x_int       = double(x_int);
        rapid_int   = double(rapid_int);
        mat_grid    = double(tensor_grid); % should be (N,M,Nt)
        
        % Need spacial dimension as first index in order to use (:) linearization
        x_int       = permute(x_int, [2 1 3]); % (M,N,Nt)
        rapid_int   = permute(rapid_int, [2 1 3]);
        
        % Enforce periodic boundary conditions
        if obj.settings.periodic_rapid 
            rapid_int   = mod(rapid_int + obj.rapid_grid(1), obj.rapid_grid(end)-obj.rapid_grid(1)) + obj.rapid_grid(1);
        end
        if obj.settings.periodic_BC 
            x_int       = mod(x_int + obj.x_grid(1), obj.x_grid(end)-obj.x_grid(1)) + obj.x_grid(1);
        end
        if obj.settings.reflective_BC     
            filter_R        = (x_int - obj.x_grid(end)) > 0; % hit right wall
            filter_L        = (x_int - obj.x_grid(1)) < 0; % hit left wall
            filter          = filter_R | filter_L;

            x_int           = (~filter) .* x_int + ...                    % non-reflected
                               filter_R.*(2*obj.x_grid(end) - x_int) + ... % reflected right
                               filter_L.*(2*obj.x_grid(1) - x_int);      % reflected left

            rapid_int        = (~filter) .* rapid_int + ...             % non-reflected
                                filter_R .* (-rapid_int) + ...         % reflected right
                                filter_L .* (-rapid_int);              % reflected right
        end
        
        
        % Get matrix representation of iFluidTensor and pemute spacial index
        % to first.
        x_g         = permute(obj.x_grid, [2 1 3]); % (M,N,Nt)
        rapid_g     = permute(obj.rapid_grid, [2 1 3]);
        
        mat_grid    = permute(mat_grid, [2 1 3]);
        mat_int     = zeros(obj.M, obj.N, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            rapid_i = rapid_int(:,:,i);
            x_i     = x_int(:,:,i);
            mat_g   = mat_grid(:,:,i);   
            
            if extrapolate
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), method);
            else
                % Set all extrapolation values to zero!
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), method, 0);
            end
           
            mat_tmp(isnan(mat_tmp)) = 0;
            mat_tmp = reshape(mat_tmp, obj.M, obj.N);
            mat_int(:,:,i) = mat_tmp;
        end
        
        % Reshape back to original indices
        mat_int = permute(mat_int, [2 1 3] );
        tensor_int = fluidcell(mat_int);
    end

end % end public methods

end % end classdef