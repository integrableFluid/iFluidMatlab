classdef LinearDiffusionSolver < iFluidSolver
    % Solves GHD Euler-equation using a first order step 
    
properties (Access = public)
    DKernel
    v_BG
    
end % end public properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = LinearDiffusionSolver(coreObj, theta_BG, Options)
        obj = obj@iFluidSolver(coreObj, Options);
        
                
        [D, v_eff] = obj.calcDiffusionKernel(theta_BG);
        obj.DKernel = D;
        obj.v_BG = v_eff;
        
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
        theta   = theta_init;
        u       = u_init;
        w       = w_init;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)
        % =================================================================
        % Purpose : Performs a single, first-order Euler step propagating
        %           the filling function theta(t) --> theta(t+dt).
        %           Note, the first order step is already implemented in 
        %           the superclass.
        % Input :   theta_prev -- Filling function at time t.
        %           u_prev     -- Position characteristic at time t.
        %           w_prev     -- Rapidity characteristic at time t.
        %           t          -- Starting time.
        %           dt         -- Length of time step.
        % Output:   theta_next -- Filling function at time t+dt.
        %           u_next     -- Position characteristic at time t+dt.
        %           w_next     -- Rapidity characteristic at time t+dt.
        % =================================================================
        [theta_next, u_next, w_next] = obj.performFirstOrderStep(theta_prev, u_prev, w_prev, t, dt);     
    end
    
    
end % end protected methods


methods (Access = public)
    
    function [D, v_eff] = calcDiffusionKernel(obj, theta_BG)
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        [rho_BG, rhoS_BG] = obj.coreObj.transform2rho(theta_BG, 0);
        
        I           = iFluidTensor(obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes);
        I.setIdentity();
        T           = -1/(2*pi)*obj.coreObj.getScatteringRapidDeriv(0, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);
        T_dr        = obj.coreObj.applyDressing(T, theta_BG, 0);
        v_eff       = obj.coreObj.calcEffectiveVelocities(theta_BG, 0, obj.x_grid, obj.rapid_grid, obj.type_grid);
        f           = obj.coreObj.getStatFactor( theta_BG );
        
        W           = rho_BG.*f.*T_dr.^2 .* abs(v_eff - transpose(v_eff));
%         w           = W*iFluidTensor(obj.rapid_w);
        w           = sum(W,1)*obj.rapid_w;
        
        D           = 1./rhoS_BG.^(2) .* (I.*w - W);
        
        D           = D.getX(1);
        v_eff       = v_eff.getX(1);
    end
    
    
    function Dk = calcDiffusiveVelocity(obj, k_array, diffFlag)
        if nargin < 3
            diffFlag = true;
        end
        
        Dk = cell(1, length(k_array));
        
        I       = iFluidTensor(obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes);
        I.setIdentity();
        
        for i = 1:length(k_array)
            Dk{i} = 1i*k_array(i)*I.*obj.v_BG + diffFlag*k_array(i)^2 * obj.DKernel/2;
        end
        
    end
    
    
    function [dtheta_t, dthetak_t, k_array] = propagateThetaTEST(obj, dtheta, t_array, diffFlag)
        % OVERLOAD
        
        %  *** First find coefficients c ***
        % Fourier transform dtheta to get g
        % Calc eigenfunctions of Dk for each relevant k --> f
        % Invert matrix of f on k and rapid (or better use \) to find c's
        
        % *** perform time evolution ***
        
        if nargin < 4
            diffFlag = true;
        end
        
        
        
        dtheta_t = cell(1, length(t_array));
        
        phi = fft( double(dtheta), [], 2 );
             
       
        P1 = phi;
        k_array = 1/(obj.x_grid(2)-obj.x_grid(1)) *(1:(obj.M))/obj.M; % add factor 0.25 for correct
        k_array = k_array - mean(k_array);
        
%         figure
%         plot(k_array, P1(obj.N/2 , :))
        
        
        DK = obj.calcDiffusiveVelocity( k_array, diffFlag);
        
        f = zeros(obj.N, obj.N, length(k_array) );
        z = zeros(1, obj.N, length(k_array) );
        c = zeros(1, obj.N, length(k_array) );
        
        for k = 1:length(k_array)
            DK_sq = permute( double(DK{k}), [1 4 2 3 5] );
            [vecs, vals] = eig(DK_sq);
            
            
            f(:,:,k) = vecs;
            z(:,:,k) = diag(vals);
            c(:,:,k) = vecs\P1(:,k); 
            
        end
        
        
        %% TEST
%         
%         gk = zeros(obj.N, 1, length(k_array));
%         
%         for i = 1:length(k_array)
%             gk(:,:,i) = f(:,:,i)*transpose(c(:,:,i));
%         end
%         
%         gk = permute( gk, [1 3 2] );
%         dt = ifft(gk, [], 2);
%         
%         figure
%         subplot(1,3,1)
%         imagesc(double(dtheta))
%         subplot(1,3,2)
%         imagesc(real(dt))
%         subplot(1,3,3)
%         imagesc(real(dt) - double(dtheta))
        

        %% Time evolution 
        
%         dtheta_t{1} = dtheta; 
        
        dthetak_t = zeros( length(t_array), obj.M );
        dthetak_t(1,:) = obj.rapid_w*sum(phi,1);
        
        for i = 1:length(t_array)
                        
            phi_t       = sum( f.*c.*exp(-z*t_array(i)) , 2 );
            phi_t       = permute( phi_t , [1 3 2] );
            
            dtheta_t{i} = iFluidTensor( ifft( phi_t , [], 2 ) );
            dthetak_t(i,:)= obj.rapid_w*sum(phi_t,1);
        end
    end
    
    
    
    
end % end public methods

end % end classdef