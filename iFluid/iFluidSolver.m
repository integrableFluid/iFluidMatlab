classdef iFluidSolver < handle
    % Superclass for solving the main GHD equation by propagating the
    % filling function, theta. The solving algorithm is abstracted leaving
    % it to the user to provide one by extending this class.
    % 
    % 
    % =========== How to extend this class ========
    %
    % ** Step 1 ** Create new solver 
    %   
    %   mySolver < iFluidSolver
    %
    %
    % ** Step 2 ** Implement the abstract methods 
    %   
    %   initialize(obj, theta_init, u_init, w_init, t_array )
    %   step(obj, theta_prev, u_prev, w_prev, t, dt)
    %
    %
    % ** Step 3 ** Write constructor calling the super-constructor
    %
    %   function obj = mySolver(coreObj, Options)
    %       obj = obj@iFluidSolver(coreObj, Options);
    %   end
    %
    %
properties (Access = protected)
    % Grid lengths
    M               = []; % number of spatial grid-points
    N               = []; % number of rapidity grid-points
    Ntypes          = []; % number of quasi-particle types
    
    % Grids (all vectors)
    rapid_grid      = [];
    x_grid          = [];
    type_grid       = [];
    rapid_w         = []; % weights for Gaussian quadrature

    coreObj         = [];
    
    % Optional parameters (default values specified here). 
    extrapFlag      = false;    % Flag for using extrapolation during propagation (useful for open systems)
    periodRapid     = false;    % Flag for periodic boundary conditions on rapidity
    calcCharac      = true;     % Flag for if characteristics are calculated
    storeNthStep    = 1;        % Store every n'th step

end % end protected properties



methods (Abstract, Access = protected)
    
    % Abstract methods must be implemented in extending class!
    
    [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array )
    [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)

end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = iFluidSolver(coreObj, Options)        
        % iFluidSolver requires an iFluidCore object for calculating
        % velocities of quasiparticles etc.
        assert( isa( coreObj, 'iFluidCore' ) )
        
        % Copy grids from iFluidCore object
        [x_grid, rapid_grid, type_grid, rapid_w] = coreObj.getGrids();

        obj.x_grid      = x_grid;
        obj.rapid_grid  = rapid_grid;
        obj.type_grid   = type_grid;
        obj.rapid_w     = rapid_w;
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = length(type_grid);
        
        obj.coreObj     = coreObj;
        
        % Copy fields of Options struct into class properties
        if ~isempty(Options)
            fn = fieldnames(Options);
            for i = 1:length(fn)                      
                if isprop(obj, fn{i}) % only copy field if defined among properties
                    eval(['obj.',fn{i},' = Options.',fn{i},';']);
                end
            end
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
        Nstored         = 1 + ceil( Nsteps/obj.storeNthStep );
        
        % Initializes cell arrays
        theta_t         = cell(1, Nstored);      
        u_t             = cell(1, Nstored);
        w_t             = cell(1, Nstored);
                                
        u_init          = fluidcell( repmat( obj.x_grid, obj.N, 1, obj.Ntypes) );
        w_init          = fluidcell( repmat( obj.rapid_grid, 1, obj.M, obj.Ntypes) );

        % Initializes the propagation, calculating and internally storing
        % any additional quantities needed for the step-function.
        [theta, u, w]   = obj.initialize(theta_init, u_init, w_init, t_array);
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
            
            if mod(n, obj.storeNthStep) == 0
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
    
    
end % end public methods


methods (Access = protected)
    
    
    function [theta_next, u_next, w_next] = performFirstOrderStep(obj, theta_prev, u_prev, w_prev, t, dt)            
        % =================================================================
        % Purpose : Performs a single, first-order Euler step propagating
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
        [v_eff, a_eff]  = obj.coreObj.calcEffectiveVelocities(theta_prev, t, obj.x_grid, obj.rapid_grid, obj.type_grid); % should be (1xNxM)
            
        x_back          = obj.x_grid - dt*v_eff;
        r_back          = obj.rapid_grid - dt*a_eff;
        
        % Use interpolation to find theta_prev at x_back, r_back and
        % assign values to theta.
        theta_next      = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
        
        if obj.calcCharac
            u_next  = obj.interpPhaseSpace(u_prev, r_back, x_back, true); % always extrapolate u
            w_next  = obj.interpPhaseSpace(w_prev, r_back, x_back, true); % always extrapolate w
        else
            u_next  = u_prev;
            w_next  = w_prev;
        end     
    end
    
    
    function tensor_int = interpPhaseSpace(obj, tensor_grid, rapid_int, x_int, extrapFlag)
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
        
        x_g         = permute(obj.x_grid, [2 1 3]); % (M,N,Nt)
        rapid_g     = permute(obj.rapid_grid, [2 1 3]);
        
        % Enforce periodic boundary conditions
        if obj.periodRapid 
            rapid_int = mod(rapid_int + obj.rapid_grid(1), obj.rapid_grid(end)-obj.rapid_grid(1)) + obj.rapid_grid(1);
        end
        
        % Get matrix representation of iFluidTensor and pemute spacial index
        % to first.
        mat_grid    = permute(mat_grid, [2 1 3]);
        mat_int     = zeros(obj.M, obj.N, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            rapid_i = rapid_int(:,:,i);
            x_i     = x_int(:,:,i);
            mat_g   = mat_grid(:,:,i);   
            
            if extrapFlag
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
  
    
end % end protected methods

end % end classdef