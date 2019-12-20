classdef iFluidSolver < handle

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

end % end protected properties



methods (Abstract, Access = protected)
    
    [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array )
    [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)

end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = iFluidSolver(coreObj, Options)        
        
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
        % Propagates filling fraction according to GHD equation.
        % Simultaneously calculates characteristic function, u, used for
        % correlation functions.
        
        Nsteps          = length(t_array) - 1;
        
        theta_t         = cell(1, Nsteps+1);
        theta_t{1}      = theta_init;
        
        u_t             = cell(1, Nsteps+1);
        u_init          = iFluidTensor( repmat( obj.x_grid, obj.N, 1, obj.Ntypes, 1) );
        u_t{1}          = u_init;

        w_t             = cell(1, Nsteps+1);
        w_init          = iFluidTensor( repmat( obj.rapid_grid, 1, 1, obj.Ntypes, 1, obj.M) );
        w_t{1}          = w_init;

        % setup for propagation
        [theta, u, w]   = obj.initialize(theta_init, u_init, w_init, t_array);

        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Time evolution progress:');
        cpb.start();   

        % Propagate theta using stepfunction
        for n = 1:Nsteps
            dt            = t_array(n+1) - t_array(n);
            [theta, u, w] = obj.step(theta, u, w, t_array(n), dt);
            theta_t{n+1}  = theta;
            u_t{n+1}      = u;
            w_t{n+1}      = w;
            
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
        % Use single Euler step
        [v_eff, a_eff]  = obj.coreObj.calcEffectiveVelocities(theta_prev, t, obj.x_grid, obj.rapid_grid, obj.type_grid); % should be (1xNxM)
            
        x_back          = obj.x_grid - dt*v_eff;
        r_back          = obj.rapid_grid - dt*a_eff;
        
        % Use interpolation to find theta_prev at x_back, r_back and
        % assign values to theta.
        theta_next      = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
        u_next          = obj.interpPhaseSpace(u_prev, r_back, x_back, true ); % always extrapolate u
        w_next          = obj.interpPhaseSpace(w_prev, r_back, x_back, true ); % always extrapolate u      
    end
    
    
    function tensor_int = interpPhaseSpace(obj, tensor_grid, rapid_int, x_int, extrapFlag)
        % This function exists because MATLAB has different syntax between
        % interp1 and interp2, and I want to choose whether i extrapolate
        % or not with a simple TRUE/FAlSE argument.
        % ASSUME function_grid is on x_grid and rapid_grid.
        % Returns function_int with same dimensions as input function_grid
        
        % rapid_int and x_int should be (N,1,Nt,1,M)
        
        % Cast to matrix form
        x_int       = double(x_int);
        rapid_int   = double(rapid_int);
        mat_grid    = double(tensor_grid); % should be (N,1,Nt,1,M)
        
        % Need spacial dimension as first index in order to use (:) linearization
        x_int       = permute(x_int, [5 1 3 4 2]); % (M,N,Nt,1,1)
        rapid_int   = permute(rapid_int, [5 1 3 4 2]);
        
        x_g         = permute(obj.x_grid, [5 1 3 4 2]); % (M,N,Nt,1,1)
        rapid_g     = permute(obj.rapid_grid, [5 1 3 4 2]);
        
        % Enforce periodic boundary conditions
        if obj.periodRapid 
            rapid_int = mod(rapid_int + obj.rapid_grid(1), obj.rapid_grid(end)-obj.rapid_grid(1)) + obj.rapid_grid(1);
        end
        
        % Get matrix representation of iFluidTensor and pemute spacial index
        % to first.
        mat_grid    = permute(mat_grid, [5 1 3 4 2]);
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
           
            mat_tmp = reshape(mat_tmp, obj.M, obj.N);
            mat_int(:,:,i) = mat_tmp;
        end
        
        % Add dummy indices and reshape back to original indices
        mat_int = permute(mat_int, [2 5 3 4 1] );
        
        tensor_int = iFluidTensor(mat_int);
    end
  
    
end % end protected methods

end % end classdef