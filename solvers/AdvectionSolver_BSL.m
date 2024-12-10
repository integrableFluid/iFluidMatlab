classdef AdvectionSolver_BSL < handle
    % Superclass for solving the main GHD equation by propagating the
    % filling function, fill. The solving algorithm is abstracted leaving
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
    %   initialize(obj, fill_init, u_init, w_init, t_array )
    %   calculateDeparturePoints(obj, fill, t, dt)
    %   update(obj, fill, x_d, r_d, v, a)
    %
    %
    % ** Step 3 ** Write constructor calling the super-constructor
    %
    %   function obj = mySolver(model, Options)
    %       obj = obj@AdvectionSolver_BSL(model, Options);
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
    SourceObj           = []; % AdvectionSource object 
    settings            = []; % struct specifying settings
    
    
end % end protected properties



methods (Abstract, Access = protected)
    
    % Abstract methods must be implemented in extending class!
    
    [fill, u, w]   = initialize(obj, fill_init, u_init, w_init, t_array)
    [x_d, r_d, v, a]= calculateDeparturePoints(obj, fill, t, dt)
    storeVelocityFields(obj, fill, x_d, r_d, v, a) % used for AM solvers

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
        addParameter(parser, 'extrapolate', false, @(x) islogical(x)); % legacy: pass as {'boundary_condition', 'open'} argument instead
        addParameter(parser, 'periodic_rapid', false, @(x) islogical(x));
        addParameter(parser, 'periodic_BC', false, @(x) islogical(x)); % legacy: pass as 'boundary_condition' argument instead
        addParameter(parser, 'reflective_BC', false, @(x) islogical(x)); % legacy: pass as 'boundary_condition' argument instead
        addParameter(parser, 'prop_characteristics', false, @(x) islogical(x));
        addParameter(parser, 'deriv_order', 1, @(x) isscalar(x) & floor(x) == x & x >=0);
        addParameter(parser, 'tol', 1e-4, @(x) x > 0);
        addParameter(parser, 'max_iter', 50, @(x) isscalar(x) & floor(x) == x & x >=0);
        addParameter(parser, 'storeNthStep', 1, @(x) isscalar(x) & floor(x) == x & x >=0);

        % parse spatial boundary conditions
        defaultBC = 'none';
        validBC = {'none','periodic','reflective','open'};
        addParameter(parser, 'boundary_conditions', defaultBC, @(x) any(validatestring(x,validBC)));


        parse(parser,varargin{:});
        
        % set settings equal to parser results
        SourceObj       = parser.Results.source;
        obj.settings    = rmfield(parser.Results, 'source');
        obj.settings    = parser.Results;

        % handle legacy boundary conditions
        if obj.settings.extrapolate
            obj.settings.boundary_conditions = 'open';
        elseif obj.settings.periodic_BC
            obj.settings.boundary_conditions = 'periodic';
        elseif obj.settings.reflective_BC
            obj.settings.boundary_conditions = 'reflective';
        end

        obj.settings    = rmfield(obj.settings, 'extrapolate');
        obj.settings    = rmfield(obj.settings, 'periodic_BC');
        obj.settings    = rmfield(obj.settings, 'reflective_BC');

        
        % handle source term (if exists)
        if ~isempty(SourceObj )
            % Inject the Solver into the Source
            SourceObj.attachSolver(obj); 

            % Store reference to Source
            obj.SourceObj = SourceObj;
            
            % cannot propagate characteristics when there is a source term
            obj.settings.prop_characteristics = false;
        end
    end
    
    
    function settings = getSettings(obj)
        settings = obj.settings;
    end

    function setBoundaryConditions(obj, boundary_conditions)
        validBC = {'none','periodic','reflective','open'};

        isValid = ismember(lower(boundary_conditions), validBC); % Case-insensitive check

        if isValid
            obj.settings.boundary_conditions = lower(boundary_conditions);
        else
            error("The value of 'boundary_conditions' is invalid. Expected input to match one of these values: %s", ...
              strjoin(validBC, ', '));
        end
    end

    
    function [fill_next, u_next, w_next] = step(obj, fill, u, w, t, dt)
        % =================================================================
        % Purpose : Propagates filling function fill one timestep dt.
        % Input :   fill   -- Filling function at time t.
        %           u       -- Lin. advection: Position characteristic
        %                      Non-lin. advection: Auxiliary variable
        %           w       -- Lin. advection: Rapidity characteristic
        %                      Non-lin. advection: Auxiliary variable
        %           t       -- Time
        %           dt      -- Time-step
        % Output:   fill_next -- Propagated filling.
        %           u_next     -- Propagated u.
        %           w_next     -- Propagated w.
        % =================================================================
        
        
        [x_d, r_d, v, a] = obj.calculateDeparturePoints(fill, t, dt);
        
        if isempty(obj.SourceObj)
            % no source term --> pure advection

            fill_next = obj.interpPhaseSpace(fill, r_d, x_d);
            obj.storeVelocityFields(fill_next, x_d, r_d, v, a);

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
            [fill_next, u_next, w_next] = obj.SourceObj.step(fill, x_d, r_d, u, w, t, dt);
            obj.storeVelocityFields(fill_next, x_d, r_d, v, a);
        end
            
    end
    
    function [fill_t, u_t, w_t] = propagateFilling(obj, fill_init, t_array)
        % =================================================================
        % Purpose : Propagates the filling function according to the GHD
        %           Euler-scale equation.
        % Input :   fill_init -- Initial filling function (fluidcell).
        %           t_array    -- Time steps for propagation
        % Output:   fill_t    -- Cell array of filling function,
        %                         with each entry corresponding to t_array.
        %           u_t        -- Cell array of position characteristics
        %                         with each entry corresponding to t_array.
        %           w_t        -- Cell array of rapidity characteristics
        %                         with each entry corresponding to t_array.
        % =================================================================
        [fill_t, u_t, w_t] = obj.propagateTheta(fill_init, t_array);
    end
    
    function [fill_t, u_t, w_t] = propagateTheta(obj, fill_init, t_array)
        % =================================================================
        % Purpose : Propagates the filling function according to the GHD
        %           Euler-scale equation.
        % Input :   fill_init -- Initial filling function (fluidcell).
        %           t_array    -- Time steps for propagation
        % Output:   fill_t    -- Cell array of filling function,
        %                         with each entry corresponding to t_array.
        %           u_t        -- Cell array of position characteristics
        %                         with each entry corresponding to t_array.
        %           w_t        -- Cell array of rapidity characteristics
        %                         with each entry corresponding to t_array.
        % =================================================================
        Nsteps          = length(t_array) - 1;
        Nstored         = 1 + ceil( Nsteps/obj.settings.storeNthStep );
        
        % Initializes cell arrays
        fill_t         = cell(1, Nstored);      
        u_t             = cell(1, Nstored);
        w_t             = cell(1, Nstored);
        
        w_init      = fluidcell( repmat( obj.rapid_grid, 1, obj.M, obj.Ntypes) );
        u_init      = fluidcell( repmat( obj.x_grid, obj.N, 1, obj.Ntypes) );
            
        % Initializes the propagation, calculating and internally storing
        % any additional quantities needed for the step-function.
        [fill, u, w]   = obj.initialize(fill_init, u_init, w_init, t_array);
        
        if ~isempty(obj.SourceObj)
            [fill, u, w] = obj.SourceObj.initialize(fill_init, u_init, w_init, t_array);
        end
        
        fill_t{1}      = fill;
        u_t{1}          = u;
        w_t{1}          = w;
        
        % Display solver settings
        if isempty(obj.SourceObj)
            source_str = 'none';
        else
            source_str = class(obj.SourceObj);
        end
        
        fprintf("Initializing BSL propagation with boundary conditions: '%s' and source: '%s'\n", ...
                obj.settings.boundary_conditions, source_str);
         
        % Initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Time evolution progress:');
        cpb.start();   

        % Propagate fill using stepfunction
        count = 2;
        for n = 1:Nsteps
            dt            = t_array(n+1) - t_array(n);
            [fill, u, w] = obj.step(fill, u, w, t_array(n), dt);
            
            if mod(n, obj.settings.storeNthStep) == 0
                fill_t{count}  = fill;
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
    

    function Vq = interpPhaseSpace(obj, V, Rq, Xq, varargin)
        % =====================================================================
        % Purpose : Interpolates the values of a fluidcell, defined on the
        %           stored rapidity and positional grids, to query points.
        %           Handles boundary conditions.
        % Input :   R, X    -- Sample rapidity (R) and position (X) grid points
        %           V       -- Sample values to interpolate (must be fluidcell 
        %                       sized [N, M, NT] )
        %           Rq, Xq  -- Query points to interpolate to (must be 
        %                       fluidcell sized [N, M, NT] )
        %           varargin-- Contains optional flags for the following
        %                       - employ extrapolation
        %                       - enforce boundary conditions
        % Output:   Vq      -- Interpolated values (returned as fluidcell)
        % =====================================================================
        
        % Define default values for the optional inputs
        defaultValues = [false, true];
    
        % Set the values of optional inputs based on varargin
        if nargin >= 5
            use_extrap = varargin{1};
        else
            use_extrap = defaultValues(1);
        end
    
        if nargin >= 6
            enforce_BC = varargin{2};
        else
            enforce_BC = defaultValues(2);
        end


        % Set interpolation method
        method = 'spline';
        if GPU_mode_on()
            method = 'cubic';
        end
        
        % Need spacial dimension as first index in order to use (:) linearization
        Xq      = permute(double(Xq), [2 1 3]);
        Rq      = permute(double(Rq), [2 1 3]);    
        V       = permute(double(V), [2 1 3]);
        Vq      = zeros(obj.M, obj.N, obj.Ntypes);


        % Enforce periodic boundary conditions of rapidity
        if obj.settings.periodic_rapid 
            Rq      = mod(Rq + obj.rapid_grid(1), obj.rapid_grid(end)-obj.rapid_grid(1)) + obj.rapid_grid(1);
        end

        % If required, enforce boundary conditions of position
        if enforce_BC && strcmpi(obj.settings.boundary_conditions, 'periodic')
            Xq      = mod(Xq + obj.x_grid(1), obj.x_grid(end)-obj.x_grid(1)) + obj.x_grid(1);

        elseif enforce_BC && strcmpi(obj.settings.boundary_conditions, 'reflective')
            filter_R= (Xq - obj.x_grid(end)) > 0;               % right boundary
            filter_L= (Xq - obj.x_grid(1)) < 0;                 % left boundary
            filter  = filter_R | filter_L;

            Xq      = (~filter) .* Xq + ...                     % non-reflected
                       filter_R.*(2*obj.x_grid(end) - Xq) + ... % reflected right
                       filter_L.*(2*obj.x_grid(1) - Xq);        % reflected left

            Rq      = (~filter) .* Rq + ...                     % non-reflected
                       filter_R .* (-Rq) + ...                  % reflected right
                       filter_L .* (-Rq);                       % reflected right

        elseif enforce_BC && strcmpi(obj.settings.boundary_conditions, 'open')
            use_extrap = true;
        end
        
        for i = 1:obj.Ntypes % interpolate for each particle type index
            Rq_i    = Rq(:,:,i);
            Xq_i    = Xq(:,:,i);
            V_i     = V(:,:,i);   
            
            if use_extrap
                % Extrapolate beyond the defined grids
                V_tmp   = interp2( obj.rapid_grid(:)', obj.x_grid(:), V_i, Rq_i(:), Xq_i(:), method);
            else
                % Set all extrapolation values to zero!
                V_tmp   = interp2( obj.rapid_grid(:)', obj.x_grid(:), V_i, Rq_i(:), Xq_i(:), method, 0);
            end
           
            V_tmp(isnan(V_tmp)) = 0;
            V_tmp       = reshape(V_tmp, obj.M, obj.N);
            Vq(:,:,i)   = V_tmp;
        end
        
        % Permute back to normal fluidcell index structure
        Vq      = permute(Vq, [2 1 3] );
        Vq      = fluidcell(Vq);
    end

        
    
    function V = calcVelocityDerivatives(obj, order, fill, t)
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration plus
        %           their temporal derivatives up to given order.
        % Input :   order   -- Max order of derivatives calculated.
        %           fill   -- Filling function (fluidcell)
        %           t       -- Time
        % Output:   v       -- cell with velocities and derivatives.
        % =================================================================
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        
        V           = cell(2, order);
        
        % calculate velocities
        Uinv        = obj.model.calcDressingOperator(fill, t);
        [veff, aeff, de_dr, dp_dr] = obj.model.calcEffectiveVelocities(fill, t, Uinv);
        
        V{1,1}      = veff;
        V{2,1}      = aeff;
        
        if order > 0  % calculate first order derivatives
            kernel      = 1/(2*pi)*obj.model.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);        

            elem        = dp_dr;
            weff        = de_dr;
            beff        = elem.*aeff;

            fill_dt    = - veff.*gradient(fill,2,obj.x_grid) - aeff.*gradient(fill,1,obj.rapid_grid);
            cell_dt     = - gradient(weff,2,obj.x_grid) - gradient(beff,1,obj.rapid_grid); 

            S_w         = -kernel*(obj.rapid_w.*fill_dt.*weff);
            S_b         = -kernel*(obj.rapid_w.*fill_dt.*beff);

            weff_dt     = obj.model.dress(S_w, Uinv);
            beff_dt     = obj.model.dress(S_b, Uinv);
            
            veff_dt     = (weff_dt - cell_dt.*veff)./elem;
            aeff_dt     = (beff_dt - cell_dt.*aeff)./elem;

            V{1,2}      = veff_dt;
            V{2,2}      = aeff_dt;
        end
           
        if order > 1 % calculate second order derivatives
            fill_dt2   = - veff_dt.*gradient(fill,2,obj.x_grid) - veff.*gradient(fill_dt,2,obj.x_grid) ...
                          - aeff_dt.*gradient(fill,1,obj.rapid_grid) - aeff.*gradient(fill_dt,1,obj.rapid_grid);

            cell_dt2    = - gradient(weff_dt,2,obj.x_grid) - gradient(beff_dt,1,obj.rapid_grid); 

            S_w2        = - kernel*(obj.rapid_w.*fill_dt2.*weff) - 2*kernel*(obj.rapid_w.*fill_dt.*weff_dt);
            S_b2        = - kernel*(obj.rapid_w.*fill_dt2.*beff) - 2*kernel*(obj.rapid_w.*fill_dt.*beff_dt);

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