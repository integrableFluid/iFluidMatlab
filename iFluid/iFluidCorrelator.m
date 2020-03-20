classdef iFluidCorrelator
    % Class for calculating Euler-scale dynamical correlation functions. 
    
    
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

    rapid_aux       = []; % auxillary grids
    type_aux        = [];
    
    TBA             = [];

    % Options
    periodRapid     = false;
    
end % end private properties


methods (Access = public)
    
    % constructor
    function obj = iFluidCorrelator(TBA, Options)       
        % iFluidCorrelator requires an iFluidCore object for TBA functions.
        assert( isa( TBA, 'iFluidCore' ) )
        
        % Copy grids from iFluidCore object
        [x_grid, rapid_grid, type_grid, rapid_w] = TBA.getGrids();

        obj.x_grid      = x_grid;
        obj.rapid_grid  = rapid_grid;
        obj.type_grid   = type_grid;
        obj.rapid_w     = rapid_w;
        
        obj.rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        obj.type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = length(type_grid);
        
        obj.TBA         = TBA;
        
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
    
    
    
    function [direct, indirect] = calc2PCorrelations( obj, theta_t, u_t, t_array, y_indices, VO1_t, VO2_0 )
        % =================================================================
        % Purpose : Calculates dynamical two-point correlation functions of
        %           the form <O1(x,t) O2(y,0)>
        % Input :   theta_t   -- cell array of filling functions. First
        %                         entry must be for time t=0.
        %           u_t       -- cell array of characteristics. First entry
        %                         must be for time t=0.
        %           t_array   -- array of times corresponding to entries in
        %                         theta_t and u_t.
        %           y_indices -- array of indices of x_grid, yielding the
        %                         desired values of y.
        %           VO1_t     -- cell array containing form factors of
        %                         operator O1 at t_array(2:end).
        %           VO2_0     -- form factor of operator O2 at time t=0.
        % Output:   direct    -- direct correlations
        %           indirect  -- indirect correlations
        % =================================================================
        assert( length(theta_t) == length(u_t) )
        assert( length(u_t) == length(t_array) )
        assert( length(t_array) > 1 )
        assert( t_array(1) == 0, "First entry of cell arrays must be at t = 0!" )
        assert( min(y_indices) >= 1, "y indices must be above one!" )
        assert( max(y_indices) <= obj.M, "y indices may not exceed length of x_grid!" )
        if length(t_array) > 2
            assert( length(VO1_t) == length(theta_t)-1 )
        end
        assert( isa( VO2_0, 'iFluidTensor' ) )
            
        
        direct      = zeros( obj.M, length(y_indices), length(t_array)-1 );
        indirect    = zeros( obj.M, length(y_indices), length(t_array)-1 );
        
        [rho_t, rhoS_t] = obj.TBA.transform2rho(theta_t, t_array);
        IM          = obj.calcStateInhomogeniety( theta_t{1}, rho_t{1} );

        
        % initialize progress bar
        Ntsteps = length(t_array) - 1;
        Nysteps = length(y_indices);
        
        fprintf('Calculating correlations for %d y-values and %d points in time.\n', Nysteps, Ntsteps )
        
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Correlation progress:');
        cpb.start();   
        
        for ti = 2:length(t_array)
            
            % Calculate gradient of characteristic u
            dudr_t = zeros(obj.N, obj.M, obj.Ntypes);
            for i = 1:obj.Ntypes
                [~, dudr_t(:,:,i)] = gradient(u_t{ti}.getType(i,'d'), 0, obj.rapid_grid);
            end
            dudr_t      = iFluidTensor( dudr_t );
        
            yc          = 1; % y_count
            
            for yi = y_indices
                                
                theta_0y    = theta_t{1}.getX(yi);
                f_0y        = obj.TBA.getStatFactor(theta_0y);
                VO2_0y      = VO2_0.getX(yi);

                corr_prod   = rhoS_t{ti}.*theta_0y.*f_0y.*VO2_0y./abs(dudr_t); 

                % Calculate source terms required for indirect propagator
                W2_temp    = rhoS_t{1}.getX(yi).*f_0y.*VO2_0y;  
                W2_temp_dr = obj.TBA.applyDressing( W2_temp, theta_0y, t_array(1) ); 
                W2_temp_sdr= W2_temp_dr - W2_temp;
                
                
                W1         = 0;
                W2         = 0;
                W3         = 0;
                
                % Iterate over x-grid, starting at x_0 where source is zero
                for xi = 1:obj.M
                    dx          = obj.x_grid(2) - obj.x_grid(1);
                    f_tx        = obj.TBA.getStatFactor(theta_t{ti}.getX(xi));
                    gamma       = obj.findRootSet( u_t{ti}.getX(xi), obj.x_grid(yi), dudr_t.getX(xi) );

                    % Calculate contribution from direct propagator
                    direct_temp = obj.interp2Gamma( corr_prod.getX(xi).*VO1_t{ti-1}.getX(xi), gamma);
                    direct(xi,yc,ti-1) = sum( sum(direct_temp,3) ,1, 'd'); % sum over gamma (rapid1 and type1)

                    
                    % ----- Calculate indirect propagator ------

                    % Update First source term, W1 
                    if isempty(gamma)
                        integr = 0;
                    else
                        Kern    = -1/(2*pi)*obj.TBA.getScatteringRapidDeriv( t_array(ti), obj.x_grid(xi), obj.rapid_grid, ...
                                                                            permute(gamma, [4 2 5 1 3]) , obj.type_grid, obj.type_aux );
                        Kern_dr = obj.TBA.applyDressing(Kern, theta_t{ti}.getX(xi), t_array(ti));
                        integr  = Kern_dr * obj.interp2Gamma(corr_prod.getX(xi) , gamma); % sum over gamma via multiplication (equal to sum over rapid2 and type2) 
                    end
                    W1          = W1 + dx*integr;
                    
                    % Calculate second source term, W2
                    W2          = - heaviside( u_t{ti}.getX(xi,'d') - obj.x_grid(yi) ) .* W2_temp_sdr;

                    % Solve for Delta (indirect propagator)
                    IM_u        = obj.interp2u( IM, u_t{ti}.getX(xi) ); % evaluate a_eff0 at x = u(t_corr, x_corr, lambda)

                    kernel      = 1/(2*pi)*obj.TBA.getScatteringRapidDeriv( t_array(ti), obj.x_grid(xi), obj.rapid_grid, ...
                                                                            obj.rapid_aux , obj.type_grid, obj.type_aux );
                    Id          = iFluidTensor( obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes, 'eye');

                    U           = Id + kernel.*transpose( obj.rapid_w.*theta_t{ti}.getX(xi) ); % CHANGED SIGN

                    vec         = rhoS_t{ti}.getX(xi) .* f_tx;
                    Ymat        = Id.*transpose(1 + 2*pi*dx*IM_u.*vec) - 2*pi*dx*IM_u.*inv(U).*transpose(vec); % vec should be transposed for consistency!!
                    
                    Delta       = Ymat\(2*pi*IM_u.*(W1+W2+W3)); % Solve integral equation through iteration
                    
                    % IMPORTANT! update W3 for next step
                    integr      = vec .* Delta;
                    integr_sdr  = obj.TBA.applyDressing( integr, theta_t{ti}.getX(xi), t_array(ti)) - integr; % should be (1xNxM)
                    W3          = W3 + dx*integr_sdr;         

                    % Calculate indirect contribution via Delta
                    indir_temp = Delta.*rho_t{ti}.getX(xi).*f_tx.*VO1_t{ti-1}.getX(xi);
                    indirect(xi,yc,ti-1) = sum( obj.rapid_w .* sum(indir_temp, 3) ,1, 'd'); % integrate over rapidity and sum over type
                end
                
                yc = yc + 1;
                
                % show progress
                Ndone = (yc-1) + (ti-2)*Nysteps;
                cpb_text = sprintf('%d/%d iterations done', Ndone, Ntsteps*Nysteps);
                cpb.setValue(Ndone/Ntsteps/Nysteps);
                cpb.setText(cpb_text);
            end
        end
        
        fprintf('\n')
    end
    
end % end public methods


methods (Access = private)

    function gamma = findRootSet(obj, u_xt, y, du_xt)
        % Finds gamma such that u(x,t,gamma) = u_xt(gamma) = y
        
        u_xt    = double(u_xt);
        du_xt   = double(du_xt);
        gamma_c = cell(1, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            % Turn u_xt into continuous, anonymous function for finding roots
            % NOTE: If fzero fails, it is often due to extrapolation. Thus, Try
            % different algorithms for interp1!
    %         u_func = @(rapid) pchip(obj.rapid_grid, u_xt, rapid);
            u_func = @(rapid) interp1(obj.rapid_grid, u_xt(:,:,i), rapid, 'linear','extrap');

            if all( du_xt(:,:,i) < 0)
                % If du_xt is monotonically decreasing with rapidity, the root
                % set contains only a single member.

                gamma_i = fzero(@(x) u_func(x) - y, 0);
            else
                % Multiple members in root set. Use sign flip to gauge how
                % many there are. 
                % NOTE: does not take into account u = y continuously
                
                % Check "bulk" of u_xt for crossings with y
                signflip    = diff( u_xt(2:end-1,:,i) - y >=0 ) ~= 0; % logical N-3 length vector indicating signflips
                
                % Check if edge-points are equal to y
                first       = abs( u_xt(1,:,i) - y ) < 1e-10;
                last        = abs( u_xt(end,:,i) - y ) < 1e-10;
                
                rapid0      = obj.rapid_grid( [first ; signflip; false; last] ); % get root-rapidities
                Nroots      = length( rapid0 ); % number of 1's in rapid0
                gamma_i     = zeros(1,Nroots);

                for j = 1:Nroots
                    gamma_i(j) = fzero(@(x) u_func(x) - y, rapid0(j));
                end    
                
                % Enforce periodic boundaries
                rapid_min = obj.rapid_grid(1);
                rapid_max = obj.rapid_grid(end);
                if obj.periodRapid 
                    gamma_i = mod(gamma_i + rapid_min, rapid_max-rapid_min) + rapid_min;
                end
                
                % Enforce uniquesnes of roots
                gamma_i     = unique(gamma_i);
            end
            
            gamma_c{i} = gamma_i;
        end
        
        % Some species might have more roots than others! To return gamma
        % as matrix, one each gamma_i must have same length. Thus, fill
        % with NaN to obtain equal length.
        maxroots    = max( cellfun(@length, gamma_c) );
        gamma       = NaN(maxroots, 1, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            Nroots      = length( gamma_c{i} );
            gamma( 1:Nroots, :, i ) = gamma_c{i};
        end
    end
    
    
    function tensor_out = interp2Gamma(obj, tensor_in, gamma)
        % gamma should be matrix of dimesions (G,M,Nt) 
        % NOTE: interp1 interpolates each column if ndims > 1
        
        if isempty(gamma)
            tensor_out = iFluidTensor(0);
            return
        end
        
        out_size    = size(tensor_in);
        out_size(1) = size(gamma,1);
        mat_out     = zeros(out_size);
        
        for i = 1:obj.Ntypes
            int = interp1(obj.rapid_grid, tensor_in.getType(i,'d'), gamma(:,:,i), 'makima');
%             int = interp1(obj.rapid_grid, tensor_in.getType(i,'d'), gamma(:,:,i), 'pchip');
%             int = interp1(obj.rapid_grid, tensor_in.getType(i,'d'), gamma(:,:,i), 'spline');
            int( isnan(int) ) = 0; % removes any NaN 
            mat_out(:,:,i) = int;
        end

        tensor_out  = iFluidTensor(mat_out);
    end
    
    
    function tensor_out = interp2u(obj, tensor_in, u_x)
        % Interpolates tensor spatially to T( lambda, u(lambda, x) ),
        % where u_x(lambda) = u(x, lambda)
        % NOTE: interp1 interpolates each column if ndims > 1

        mat_out = zeros(obj.N, 1, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            tens_i  = tensor_in.getType(i,'d');
            ux_i    = u_x.getType(i,'d');
            
            % make sure spatial index is first
            tens_i  = permute(tens_i, [2 1]);
            x_g     = permute(obj.x_grid, [2 1]);
            
%             mat_int = interp1(x_g, tens_i, ux_i, 'spline'); % (N,N)
            mat_int = interp1(x_g, tens_i, ux_i, 'makima'); % (N,N)
            
            % Since u(x, lambda) must be evaluated at same rapidities as
            % the tensor T, one has to take the diagonal.            
            mat_out(:,:,i) = diag( mat_int );
        end

        tensor_out = iFluidTensor(mat_out);
    end
    
    
    function IM = calcStateInhomogeniety( obj, theta_0, rho_0 )
        % Calculate the inhomogeneity of the initial state by taking a
        % numerical derivative of theta_0
        
        x_g     = permute(obj.x_grid, [2 1]);
        x_fine  = linspace( x_g(1), x_g(end), 3*(length(x_g)-1)+1 );
        
        IM      = zeros(size(theta_0));
        
        for i = 1:obj.Ntypes
            t0i     = theta_0.getType(i,'d');
            r0i     = rho_0.getType(i,'d');
            
            % make sure spatial index is first
            t0i     = permute(t0i, [2 1]);
            r0i     = permute(r0i, [2 1]);
            
            % interpolate to a finer grid to help improve gradient accuracy
            t_int   = interp1(x_g, t0i, x_fine, 'makima');
            r_int   = interp1(x_g, r0i, x_fine, 'makima');
            
            [~,grad]= gradient(t_int, 0, x_fine);
            IM_temp = grad./(eps + 2*pi*r_int.*obj.TBA.getStatFactor(t_int));
            IM_temp(isnan(IM_temp)) = 0;
            IM_temp = permute(IM_temp, [2 1]);
            
            IM(:,:,i) = IM_temp(:, 1:3:end);
        end
 
        IM = iFluidTensor(IM);
    end
    
end % end private methods
    
end % end classdef 
