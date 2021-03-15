classdef LargeDeviationModule
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
    function obj = LargeDeviationModule(TBA, Options)       
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
    
    
    function [direct, indirect] = calc2PCorrelations( obj, theta_t, U_t, VO_t, t_array, y_indices)
        % =================================================================
        % Purpose : Calculates dynamical two-point correlation functions of
        %           the form <O1(x,t1) O2(y,t2)>
        % Input :   theta_t   -- (2,Nt)-sized cell array of filling 
        %                         functions.
        %                         First row is theta at times t1,
        %                         Second row is theta at times t2.
        %           U_t       -- (1,Nt)-sized cell array of characteristic 
        %                         functions. The i'th entry should be 
        %                         U(t_array(1,i), t_array(2,i))
        %           t_array   -- (2,Nt)-sized array of times.
        %                         First row contains times t1.
        %                         Second row contains times t2.
        %           VO_t      -- (2,Nt)-sized cell array of form factors.
        %                         First row is V1 at times t1,
        %                         Second row is U at times t2.
        %           y_indices -- array of indices of x_grid, yielding the
        %                         desired values of y.
        % Output:   direct    -- direct correlations
        %           indirect  -- indirect correlations
        % =================================================================           
        
        
        Nt          = size(theta_t, 2);
        
        direct      = zeros(obj.M, length(y_indices), Nt);
        indirect    = zeros(obj.M, length(y_indices), Nt);
        
        
        [rho_t(1,:), rhoS_t(1,:)] = obj.TBA.transform2rho(theta_t(1,:), t_array(1,:));
        [rho_t(2,:), rhoS_t(2,:)] = obj.TBA.transform2rho(theta_t(2,:), t_array(2,:));
        
        IM          = obj.calcStateInhomogeniety( theta_t(2,:), rho_t(2,:) );

        
        % initialize progress bar
        Ntsteps = Nt;
        Nysteps = length(y_indices);
        
        fprintf('Calculating correlations for %d y-values and %d points in time.\n', Nysteps, Ntsteps )
        
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Correlation progress:');
        cpb.start();   
        
        for ti = 1:Nt
            
            % Calculate gradient of characteristic U
            dUdr = zeros(obj.N, obj.M, obj.Ntypes);
            for i = 1:obj.Ntypes
                [~, dUdr(:,:,i)] = gradient(U_t{ti}.getType(i,'d'), 0, obj.rapid_grid);
            end
            dUdr = iFluidTensor( dUdr );
        
            yc   = 1; % y_count
            
            for yi = y_indices
                % Scan over x to calculate correlations
                [dir, indir] = obj.integrateLeftToRight(theta_t{1,ti}, theta_t{2,ti}, ...
                                                        rhoS_t{1,ti}, rhoS_t{2,ti}, ...
                                                        VO_t{1,ti}, VO_t{2,ti}, ...
                                                        U_t{ti}, dUdr, ...
                                                        t_array(1,ti), t_array(2,ti), ...
                                                        IM{ti}, yi);            
  
                direct(:,yc,ti) = dir;
                indirect(:,yc,ti) = indir;
                
                yc = yc + 1;
                
                % show progress
                Ndone = (yc-1) + (ti-1)*Nysteps;
                cpb_text = sprintf('%d/%d iterations done', Ndone, Ntsteps*Nysteps);
                cpb.setValue(Ndone/Ntsteps/Nysteps);
                cpb.setText(cpb_text);
            end
        end
        
        fprintf('\n')
    end
    
    
    function [U, W] = calcCharacteristics(obj, theta_t, t_array, t_0_idx)
        % =================================================================
        % Purpose : Calculates the characteristics with starting points 
        %           t_array(t_0_idx).
        % Input :   theta_t    -- Filling function at times t in t_array.
        %           t_array    -- Array of times t (gives propagation step).
        %           t_0_idx    -- Starting indices of t_array.
        % Output:   U          -- Position characteristics.
        %           W          -- Rapidity characteristics.
        % =================================================================
        
        % Initialize characteristics
        U_init      = iFluidTensor( repmat( obj.x_grid, obj.N, 1, obj.Ntypes) );
        W_init      = iFluidTensor( repmat( obj.rapid_grid, 1, obj.M, obj.Ntypes) );
        
        U           = cell(length(t_0_idx), length(t_array));
        W           = cell(length(t_0_idx), length(t_array));
        
        for i = 1:length(t_0_idx)
            U{i,t_0_idx(i)} = U_init;
            W{i,t_0_idx(i)} = W_init;
        end
        
        % Propagate characteristic
        for n = 1:length(t_array)-1
            dt          = t_array(n+1) - t_array(n);
            t_cur       = t_array(n);
            theta_mid   = (theta_t{n+1} + theta_t{n})/2;
            
            [v_e, a_e]  = obj.TBA.calcEffectiveVelocities(theta_mid, t_cur+dt/2);
            x_mid       = obj.x_grid - 0.5*dt*v_e;
            r_mid       = obj.rapid_grid - 0.5*dt*a_e; 

            % Interpolate v_eff and a_eff to midpoint coordinates x' and rapid' 
            v_mid       = obj.interpPhaseSpace( v_e, r_mid, x_mid, true );
            a_mid       = obj.interpPhaseSpace( a_e, r_mid, x_mid, true );
            
            % From midpoint velocity calculate fullstep x and rapid translation
            x_back      = obj.x_grid - dt*v_mid;
            r_back      = obj.rapid_grid - dt*a_mid;
            
            % propagate all characterstics
            for i = 1:length(t_0_idx)
                U_prev = U{i,n};
                W_prev = W{i,n};
                
                if ~isempty(U_prev) && ~isempty(W_prev)
                    U{i,n+1} = obj.interpPhaseSpace(U_prev, r_back, x_back, true);
                    W{i,n+1} = obj.interpPhaseSpace(W_prev, r_back, x_back, true); 
                end
            end          
            
        end    
    end
    
    
    function epsilon = solveFlowEquation(obj)
        % Work in progress!! A couple of notes regarding this function:
        %
        % - For a better notation, we use T instead of t'
        % - tau_array contains tau in [0 t]
        % - For now, assume T = t. In order to get the characteristic from
        %   tau to T, we need to know all the states inbetween, so we might
        %   as well compute all the epsilons inbetween while were all it.
        % - Instead of incorporating s into the existing iFluidTensor
        %   structure, we treat it like a time index, namely by storing 
        %   iFluidTensor objects in a cell array. From now on, s will be
        %   the first index of the cell array, while time will be the
        %   second index. Once this function is working, I might go back
        %   and polish the implementation a bit.
        
        ds_array    = diff(s_array); 
            
        error_rel   = 1;
        count       = 0;
        eps_old     = cell(1:length(ds_array)); % dont update eps0
        eps_old(:)  = {eps0}; % set starting guess as for all s>0 as eps(s=0) 

        while any(error_rel > obj.tolerance) & count < obj.maxcount

            eps = eps_0 + ds_array.*cumsum( flowIteration(eps_old) );


            % calculate error
            v1          = flatten(eps);
            v2          = flatten(eps_old);

            sumeff      = sum( v1.^2 ,1);            
            error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
            eps_old     = eps;

            count       = count+1;
        end
        
        
        epsilon = [eps0, eps];
        
        
        
        function eps_next = flowIteration(eps_prev)
            % I put all of this in a nested function to avoid clutter in
            % the main iteration loop. Will need to test if this negatively
            % affects performance. Possibly there is a large overhead, but
            % it remains to be seen ....
            
            
            % When evaluating the couplings, should one take the starting
            % or the ending time?? I assume the ending time here ...
            
            % First index is start time (t1), second index is end time (t2)
            % Note: t2 >= t1
            Nt          = length(t_array);
            
            theta_tt    = cell(Nt, Nt);
            rho_tt      = cell(Nt, Nt);
            rhoS_tt     = cell(Nt, Nt);
            U_tt        = cell(Nt, Nt);
            eps_next    = cell(Nt, Nt);
            
            
            % Calculate theta and rho for all time combinations
            for i = 1:Nt % starting times
                for j = i:Nt % ending times
                    theta           = obj.TBA.calcFillingFraction(eps_prev{i,j});
                    [rho, rhoS]     = obj.TBA.transform2rho(theta, t_array(j));
                    
                    theta_tt{i,j}   = theta;
                    rho_tt{i,j}     = rho;
                    rhoS_tt{i,j}    = rhoS;
                    U_tt{i,j}       = U;
                end
                
                
                U_tt(i,:) = obj.calcCharacteristics(theta_tt(i,i:Nt), t_array(i:Nt), 1);
            end
            
            TERM2       = 0;
            
            for i = 1:Nt
                for j = i:Nt
                    % Get characteristic and calculate derivative
                    U       = U_tt{i, j}; 
                    dUdr    = zeros(obj.N, obj.M, obj.Ntypes);
                    for k = 1:obj.Ntypes
                        [~, dUdr(:,:,k)] = gradient(U.getType(k,'d'), 0, obj.rapid_grid);
                    end
                    dUdr    = iFluidTensor( dUdr );


                    % Calculate the quantities needed for the first term
                    TERM1 = 0;
                    tau_star = obj.findRootSetTime(U, obj.x_grid(y_idx), t_array(i)); % output has size (Nroots,1,Ntypes)

                    % NOTE!!: For now we completely ignore multiple types in
                    % tau_star. Let's just get it to work for hard rods and LL
                    % first.

                    for j = 1:size(tau_star,1) 
                        % interpolate theta to the time tau_star
                        theta_star  = obj.interpTime(tau_array, theta_tt.getX(y_idx), tau_star(j,:,:));

                        v_eff       = obj.TBA.calcEffectiveVelocities(theta_star, tau_star(j,:,:));
                        h           = obj.TBA.getOneParticleEV( c_idx, tau_star(j,:,:), obj.x_grid(y_idx), obj.rapid_grid);               
                        hn_dr       = obj.TBA.applyDressing(h, theta_star, tau_star(j,:,:));


                        % How to take sgn() of v_eff?? Cant take eigenvalue, so
                        % I guess its just the regular sign function (????)
                        TERM1 = TERM1 - sign(double(v_eff)).*hn_dr;
                    end




                    % Calculating the remaining quantities for the second term
                    [~,~,~,Vj] = obj.TBA.calcCharges(c_idx, theta_tt{i}, tau_array(i), true);
                    IM = obj.calcStateInhomogeniety( theta_tt{i}, rho_tau{i} );


                    [~,~,Delta] = obj.integrateLeftToRight(theta_T, theta_tt{i}, ...
                                                           rhoS_T, rhoS_tau{i}, ...
                                                           Vj, Vj, ...
                                                           U, dUdr, ...
                                                           T, tau_array(i), ...
                                                           IM, y_idx);


                    TERM2 = TERM2 - dt*Delta; % one more step in integration
                    eps_next{i} = TERM1 + TERM2;
                end
            end
            
        end
        
        
    end
    
    
    function tensor_int = interpTime(obj, t_array, tensor_t, t_int)
        % public for the purpose of testing. Make private later
        % First version. Very naiive. Assume t_int is inside t_array
        % endpoints. Then find two entries in t_array around t_int and do a
        % linear interpolation between tensor at these two points.
        
        assert( all(t_int > t_array(1)) )
        assert( all(t_int < t_array(end)) )
        
        for i = 1:length(t_int)
            % Find entry in t_array lower than t_int
            [t_low, idx]    = max( t_array(t_array < t_int(i)) );
            t_high          = t_array(idx+1);
            
            dt              = t_high - t_low; % difference in grid
            s               = (t_high - t_int(i))/dt; % fractional difference to high
            
            tensor_temp     = (1-s)*tensor_t{idx+1} + s*tensor_t{idx};
            
            if length(t_int) == 1
                tensor_int = tensor_temp;
            else
                tensor_int{i} = tensor_temp;
            end
            
        end
        
    end
    
     
end % end public methods


methods (Access = private)

    function gamma = findRootSetRapid(obj, u_xt, y, du_xt)
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
    
    
    function tau = findRootSetTime(obj, U, y, tstart_array)
        % TODO: TEST THIS FUNCTION!
        
        % Finds tau such that U(x,t,lambda,tau) = y
        % Input U should be cell-array with all U having same ending time
        % but variable starting time.
        
        tau_cell = cell(1, obj.Ntypes);

        % Unpack U (cell-array) into 4d-array (starttime,rapid,space,type)
        U_arr = zeros(obj.N, obj.M, obj.Ntypes, length(U));
        for i = 1:length(U)
            U_arr = double(U{i});
        end
        U_arr = permute(U_arr, [4, 1, 2, 3]);
        [~,dU_dt] = gradient(U_arr, tstart_array);
        
        for i = 1:obj.Ntypes
            U_i = U_arr(:,:,:,i); % get U for this type
            dU_i = dU_dt(:,:,:,i); % get dU_dt for this type
            
            % Turn U into continuous, anonymous function for finding roots
            % NOTE: If fzero fails, it is often due to extrapolation. Thus, Try
            % different algorithms for interp1!
            U_func = @(t) interp1(tstart_array, U_i, t, 'linear','extrap');

            if all( dU_i < 0) || all( dU_i > 0)
                % If du_dt is monotonic, the root set contains only a 
                % single member.

                tau_i = fzero(@(x) U_func(x) - y, 0);
            else
                % Multiple members in root set. Use sign flip to gauge how
                % many there are. 
                % NOTE: does not take into account U = y continuously
                
                % Check "bulk" of u for crossings with y
                signflip    = diff( U_i(2:end-1,:,:) - y >=0 ) ~= 0; % logical N-3 length vector indicating signflips
                
                % Check if edge-points are equal to y
                first       = abs( U_i(1,:,:) - y ) < 1e-10;
                last        = abs( U_i(end,:,:) - y ) < 1e-10;
                
                t_flip      = tstart_array( [first ; signflip; false; last] ); 
                Nroots      = length( t_flip ); % number of 1's in t_flip
                tau_i       = zeros(1,Nroots);

                for j = 1:Nroots
                    tau_i(j) = fzero(@(x) U_func(x) - y, t_flip(j));
                end    
                
                % Enforce uniquesnes of roots
                tau_i     = unique(tau_i);
            end
            
            tau_cell{i} = tau_i;
        end
        
        % Some species might have more roots than others! To return gamma
        % as matrix, one each gamma_i must have same length. Thus, fill
        % with NaN to obtain equal length.
        maxroots    = max( cellfun(@length, tau_cell) );
        tau         = NaN(maxroots, 1, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            Nroots      = length( tau_cell{i} );
            tau( 1:Nroots, :, i ) = tau_cell{i};
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
    
    
    function IM_t = calcStateInhomogeniety( obj, theta_t, rho_t )
        % Calculate the inhomogeneity of the initial state by taking a
        % numerical derivative of theta
        
        x_g     = permute(obj.x_grid, [2 1]);
        x_fine  = linspace( x_g(1), x_g(end), 3*(length(x_g)-1)+1 );
        
        IM_t    = cell(1, length(theta_t));
        
        for j = 1:length(theta_t)    
            theta = theta_t{j};
            rho = rho_t{j};
            IM  = zeros(size(theta));

            for i = 1:obj.Ntypes
                t0i     = theta.getType(i,'d');
                r0i     = rho.getType(i,'d');

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

            IM_t{j} = iFluidTensor(IM);
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
        
        tensor_int = iFluidTensor(mat_int);
    end
    
    
    function [direct, indirect, Delta] = integrateLeftToRight(obj, ...
                                                       theta_t1, theta_t2, ...
                                                       rhoS_t1, rhoS_t2, ...
                                                       VO_t1, VO_t2, ...
                                                       U, dUdr, ...
                                                       t1, t2, ...
                                                       IM, yi)
        % Calculate <O1(x,t1) O2(y,t2)> for all values of x by integrating
        % from left to right.
        
        % Calculate quantities constant for all values of x
        rho_t1      = theta_t1.*rhoS_t1;
        f_t2_y      = obj.TBA.getStatFactor(theta_t2.getX(yi));
        corr_prod   = rhoS_t1.*theta_t2.getX(yi).*f_t2_y.*VO_t2.getX(yi)./abs(dUdr);
        

        % Calculate source terms required for indirect propagator
        W2_temp    = rhoS_t2.getX(yi).*f_t2_y.*VO_t2.getX(yi);  
        W2_temp_dr = obj.TBA.applyDressing( W2_temp, theta_t2.getX(yi), t2 ); 
        W2_temp_sdr= W2_temp_dr - W2_temp;
        
        
        W1      = 0;
        W2      = 0;
        W3      = 0;
        
        direct  = zeros(obj.M, 1);
        indirect= zeros(obj.M, 1);
        Delta   = zeros(obj.N, obj.M, obj.Ntypes);

        for xi = 1:obj.M
            dx          = obj.x_grid(2) - obj.x_grid(1);
            f_t1_x      = obj.TBA.getStatFactor(theta_t1.getX(xi));
            gamma       = obj.findRootSetRapid( U.getX(xi), obj.x_grid(yi), dUdr.getX(xi) );

            % ----- Calculate direct propagator ------
            direct_temp = obj.interp2Gamma( corr_prod.getX(xi).*VO_t1.getX(xi), gamma);
            direct(xi) = sum( sum(direct_temp,3) ,1, 'd'); % sum over gamma (rapid1 and type1)


            % ----- Calculate indirect propagator ------

            % Update First source term, W1 
            if isempty(gamma)
                integr = 0;
            else
                Kern    = -1/(2*pi)*obj.TBA.getScatteringRapidDeriv( t1, obj.x_grid(xi), obj.rapid_grid, ...
                                                                    permute(gamma, [4 2 5 1 3]) , obj.type_grid, obj.type_aux );
                Kern_dr = obj.TBA.applyDressing(Kern, theta_t1.getX(xi), t1);
                integr  = Kern_dr * obj.interp2Gamma(corr_prod.getX(xi) , gamma); % sum over gamma via multiplication (equal to sum over rapid2 and type2) 
            end
            W1          = W1 + dx*integr;

            % Calculate second source term, W2
            W2          = - heaviside( U.getX(xi,'d') - obj.x_grid(yi) ) .* W2_temp_sdr;

            % Solve for Delta (indirect propagator)
            IM_u        = obj.interp2u( IM, U.getX(xi) ); % evaluate a_eff0 at x = u(t_corr, x_corr, lambda)

            kernel      = 1/(2*pi)*obj.TBA.getScatteringRapidDeriv( t1, obj.x_grid(xi), obj.rapid_grid, ...
                                                                    obj.rapid_aux , obj.type_grid, obj.type_aux );
            Id          = iFluidTensor( obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes, 'eye');

            Umat        = Id + kernel.*transpose( obj.rapid_w.*theta_t1.getX(xi) ); 

            vec         = rhoS_t1.getX(xi).*f_t1_x;
            Ymat        = Id.*(1 + 2*pi*dx*IM_u.*vec) - 2*pi*dx*IM_u.*inv(Umat).*transpose(vec); % vec should be transposed for consistency!!

            Delta_x     = Ymat\(2*pi*IM_u.*(W1+W2+W3)); % Solve integral equation through iteration
            Delta(:,xi,:)= double(Delta_x);
            
            % IMPORTANT! update W3 for next step
            integr      = vec .* Delta_x;
            integr_sdr  = obj.TBA.applyDressing( integr, theta_t1.getX(xi), t1) - integr; % should be (1xNxM)
            W3          = W3 + dx*integr_sdr;         

            % Calculate indirect contribution via Delta
            indir_temp  = Delta_x.*rho_t1.getX(xi).*f_t1_x.*VO_t1.getX(xi);
            indirect(xi)= sum( obj.rapid_w .* sum(indir_temp, 3) ,1, 'd'); % integrate over rapidity and sum over type
        end
        
    end
    
end % end private methods
    
end % end classdef 
