classdef AdvectionSolver_BSL < handle
    % Superclass for solving the main GHD equation by propagating the
    % filling function, theta. The solving algorithm is abstracted leaving
    % it to the user to provide one by extending this class.
    % 
    % 
    % =========== How to extend this class ========
    %
    % ** Step 1 ** Create new solver 
    %   
    %   mySolver < AdvectionSolver_BSL
    %
    %
    % ** Step 2 ** Implement the abstract methods 
    %   
    %   initialize(obj, theta_init, u_init, w_init, t_array )
    %   calculateDeparturePoints(obj, theta, t, dt)
    %   update(obj, theta, x_d, r_d, v, a)
    %
    %
    % ** Step 3 ** Write constructor calling the super-constructor
    %
    %   function obj = mySolver(model, Options)
    %       obj = obj@iFluidSolver(model, Options);
    %   end
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
    source              = []; % AdvectionSource object 
    settings            = []; % struct specifying settings
    
    
end % end protected properties



methods (Abstract, Access = protected)
    
    % Abstract methods must be implemented in extending class!
    
    [theta, u, w]   = initialize(obj, theta_init, u_init, w_init, t_array)
    [x_d, r_d, v, a]= calculateDeparturePoints(obj, theta, t, dt)
    storeVelocityFields(obj, theta, x_d, r_d, v, a) % used for AM solvers

end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = AdvectionSolver_BSL(model, varargin)        
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

        addParameter(parser, 'source', [], @(x) isa( x, 'AdvectionSource' ));
        addParameter(parser, 'extrapolate', false, @(x) islogical(x));
        addParameter(parser, 'periodic_rapid', false, @(x) islogical(x));
        addParameter(parser, 'periodic_BC', false, @(x) islogical(x));
        addParameter(parser, 'reflective_BC', false, @(x) islogical(x));
        addParameter(parser, 'prop_characteristics', false, @(x) islogical(x));
        addParameter(parser, 'deriv_order', 1, @(x) isscalar(x) & floor(x) == x & x >=0);
        addParameter(parser, 'tol', 1e-4, @(x) x > 0);
        addParameter(parser, 'max_iter', 50, @(x) isscalar(x) & floor(x) == x & x >=0);
        addParameter(parser, 'storeNthStep', 1, @(x) isscalar(x) & floor(x) == x & x >=0);
        

        parse(parser,varargin{:});
        
        % set settings equal to parser results
        obj.source      = parser.Results.source;
        obj.settings    = rmfield(parser.Results, 'source');
        obj.settings    = parser.Results;
        
        
        if ~isempty(obj.source)
            % copy settings from solver to source term
            obj.source.copySettings(obj.settings);
            
            % cannot propagate characteristics when there is a source term
            obj.settings.prop_characteristics = false;
        end
    end
    
    
    function settings = getSettings(obj)
        settings = obj.settings;
    end

    
    function [theta_next, u_next, w_next] = step(obj, theta, u, w, t, dt)
        % =================================================================
        % Purpose : Propagates filling function theta one timestep dt.
        % Input :   theta   -- Filling function at time t.
        %           u       -- Lin. advection: Position characteristic
        %                      Non-lin. advection: Auxiliary variable
        %           w       -- Lin. advection: Rapidity characteristic
        %                      Non-lin. advection: Auxiliary variable
        %           t       -- Time
        %           dt      -- Time-step
        % Output:   theta_next -- Propagated filling.
        %           u_next     -- Propagated u.
        %           w_next     -- Propagated w.
        % =================================================================
        
        
        [x_d, r_d, v, a] = obj.calculateDeparturePoints(theta, t, dt);
        
        if isempty(obj.source)
            % no source term --> pure advection

            theta_next = obj.interpPhaseSpace(theta, r_d, x_d, obj.settings.extrapolate);
            obj.storeVelocityFields(theta_next, x_d, r_d, v, a);

            if obj.settings.prop_characteristics
                % return characteristics X(t+dt, 0)
                u_next = obj.interpPhaseSpace(u, r_n, x_n, true);
                w_next = obj.interpPhaseSpace(w, r_n, x_n, true);
            else
                % return characteristics X(t+dt, t)
                u_next = x_d;
                w_next = r_d;
            end    
        
        else
            % source term --> non-linear advection step
            % u_next and w_next are auxiliary variables
            [theta_next, u_next, w_next] = obj.source.step(theta, x_d, r_d, u, w, t, dt);
            obj.storeVelocityFields(theta_next, x_d, r_d, v, a);
        end
            
    end
    
    
    function [theta_t, u_t, w_t] = propagateTheta(obj, theta_init, t_array)
        % =================================================================
        % Purpose : Propagates the filling function according to the GHD
        %           Euler-scale equation.
        % Input :   theta_init -- Initial filling function (iFluidTensor).
        %           t_array    -- Time steps for propagation
        % Output:   theta_t    -- Cell array of filling function,
        %                         with each entry corresponding to t_array.
        %           u_t        -- Cell array of position characteristics
        %                         with each entry corresponding to t_array.
        %           w_t        -- Cell array of rapidity characteristics
        %                         with each entry corresponding to t_array.
        % =================================================================
        Nsteps          = length(t_array) - 1;
        Nstored         = 1 + ceil( Nsteps/obj.settings.storeNthStep );
        
        % Initializes cell arrays
        theta_t         = cell(1, Nstored);      
        u_t             = cell(1, Nstored);
        w_t             = cell(1, Nstored);
        
        w_init      = fluidcell( repmat( obj.rapid_grid, 1, obj.M, obj.Ntypes) );
        u_init      = fluidcell( repmat( obj.x_grid, obj.N, 1, obj.Ntypes) );
            
        % Initializes the propagation, calculating and internally storing
        % any additional quantities needed for the step-function.
        [theta, u, w]   = obj.initialize(theta_init, u_init, w_init, t_array);
        
        if ~isempty(obj.source)
            [theta, u, w] = obj.source.initialize(theta_init, u_init, w_init, t_array);
        end
        
        theta_t{1}      = theta;
        u_t{1}          = u;
        w_t{1}          = w;
        
         
        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Time evolution progress:');
        cpb.start();   

        % Propagate theta using stepfunction
        count = 2;
        for n = 1:Nsteps
            dt            = t_array(n+1) - t_array(n);
            [theta, u, w] = obj.step(theta, u, w, t_array(n), dt);
            
            if mod(n, obj.settings.storeNthStep) == 0
                theta_t{count}  = theta;
                u_t{count}      = u;
                w_t{count}      = w;
                
                count           = count + 1;
            end
            
            % show progress
            cpb_text = sprintf('%d/%d steps evolved', n, Nsteps);
            cpb.setValue(n/Nsteps);
            cpb.setText(cpb_text);
        end
        fprintf('\n')
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
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), 'spline');
            else
                % Set all extrapolation values to zero!
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), 'spline', 0);
            end
           
            mat_tmp(isnan(mat_tmp)) = 0;
            mat_tmp = reshape(mat_tmp, obj.M, obj.N);
            mat_int(:,:,i) = mat_tmp;
        end
        
        % Reshape back to original indices
        mat_int = permute(mat_int, [2 1 3] );
        tensor_int = fluidcell(mat_int);
    end
        
    
    function V = calcVelocityDerivatives(obj, order, theta, t)
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration plus
        %           their temporal derivatives up to given order.
        % Input :   order   -- Max order of derivatives calculated.
        %           theta   -- Filling function (fluidcell)
        %           t       -- Time
        % Output:   v       -- cell with velocities and derivatives.
        % =================================================================
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        
        V           = cell(2, order);
        
        % calculate velocities
        Uinv        = obj.model.calcDressingOperator(theta, t);
        [veff, aeff, de_dr, dp_dr] = obj.model.calcEffectiveVelocities(theta, t, Uinv);
        
        V{1,1}      = veff;
        V{2,1}      = aeff;
        
        if order > 0  % calculate first order derivatives
            kernel      = 1/(2*pi)*obj.model.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);        

            elem        = dp_dr;
            weff        = de_dr;
            beff        = elem.*aeff;

            theta_dt    = - veff.*gradient(theta,2,obj.x_grid) - aeff.*gradient(theta,1,obj.rapid_grid);
            cell_dt     = - gradient(weff,2,obj.x_grid) - gradient(beff,1,obj.rapid_grid); 

            S_w         = -kernel*(obj.rapid_w.*theta_dt.*weff);
            S_b         = -kernel*(obj.rapid_w.*theta_dt.*beff);

            weff_dt     = obj.model.dress(S_w, Uinv);
            beff_dt     = obj.model.dress(S_b, Uinv);
            
            veff_dt     = (weff_dt - cell_dt.*veff)./elem;
            aeff_dt     = (beff_dt - cell_dt.*aeff)./elem;

            V{1,2}      = veff_dt;
            V{2,2}      = aeff_dt;
        end
           
        if order > 1 % calculate second order derivatives
            theta_dt2   = - veff_dt.*gradient(theta,2,obj.x_grid) - veff.*gradient(theta_dt,2,obj.x_grid) ...
                          - aeff_dt.*gradient(theta,1,obj.rapid_grid) - aeff.*gradient(theta_dt,1,obj.rapid_grid);

            cell_dt2    = - gradient(weff_dt,2,obj.x_grid) - gradient(beff_dt,1,obj.rapid_grid); 

            S_w2        = - kernel*(obj.rapid_w.*theta_dt2.*weff) - 2*kernel*(obj.rapid_w.*theta_dt.*weff_dt);
            S_b2        = - kernel*(obj.rapid_w.*theta_dt2.*beff) - 2*kernel*(obj.rapid_w.*theta_dt.*beff_dt);

            weff_dt2    = obj.model.dress(S_w2, Uinv);
            beff_dt2    = obj.model.dress(S_b2, Uinv);

            veff_dt2    = (weff_dt2 - 2*veff_dt.*cell_dt - veff.*cell_dt2)./elem;
            aeff_dt2    = (beff_dt2 - 2*aeff_dt.*cell_dt - aeff.*cell_dt2)./elem; 

            V{1,3}      = veff_dt2;
            V{2,3}      = aeff_dt2;
        end
        
        if order > 2
            error('Higher order derivatives not yet implemented!')
        end
        
    end
        

end % end public methods

end % end classdef