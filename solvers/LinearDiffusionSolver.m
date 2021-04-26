classdef LinearDiffusionSolver < iFluidSolver
    % Solves GHD propagation eq. with linearized diffusion for a
    % perturbation on top of a homogeneous background state.
    
properties (Access = public)
    DKernel = []; % diffusion kernel
    v_BG    = []; % effective velocity from background
    
end % end public properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = LinearDiffusionSolver(coreObj, theta_BG, Options)
        obj         = obj@iFluidSolver(coreObj, Options);
                
        [D, v_eff]  = obj.calcDiffusionKernel(theta_BG);
        obj.DKernel = D;
        obj.v_BG    = v_eff;
    end
    
    
end % end public methods


methods (Access = protected)

    % Implementation of abstract functions
    
    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        theta = theta_init;
        u = u_init;
        w = w_init;
    end
      
    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)
        % Not used
    end
    
    
end % end protected methods


methods (Access = public)
    
    function [D, v_eff] = calcDiffusionKernel(obj, theta_BG)
        % =================================================================
        % Purpose : Calculate diffusion kernel from background state.
        % Input :   theta_BG -- Filling function of background state.
        % Output:   D        -- Linearized diffusion kernel
        %           v_eff    -- Effective velocity from background
        % =================================================================
        
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        [rho_BG, rhoS_BG] = obj.coreObj.transform2rho(theta_BG, 0);
        
        % Prepare required quantities
        I           = fluidcell.eye(obj.N, obj.Ntypes)./obj.rapid_w;
        T           = -1/(2*pi)*obj.coreObj.getScatteringRapidDeriv(0, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);
        T_dr        = obj.coreObj.applyDressing(T, theta_BG, 0);
        v_eff       = obj.coreObj.calcEffectiveVelocities(theta_BG, 0, obj.x_grid, obj.rapid_grid, obj.type_grid);
        f           = obj.coreObj.getStatFactor( theta_BG );
        
        % Calculate diffusion kernel
        W           = rho_BG.*f.*T_dr.^2 .* abs(v_eff - v_eff.t());
        w           = sum( transpose(W.*obj.rapid_w) , 4);
        D           = rhoS_BG.^(-2).*(I.*w - W).*sqrt(obj.rapid_w.*permute(obj.rapid_w, [4 2 3 1]));
        
    end
    
    
    function Dk = calcDiffusiveVelocity(obj, k_array, diffFlag)
        % =================================================================
        % Purpose : Calculate diffusion operator for each mode.
        % Input :   k_array -- Array of momentum modes.
        %           diffFlag-- if 0, disable diffusion.
        % Output:   Dk      -- Cell array of diff. ops. for each mode k.
        % =================================================================
        
        if nargin < 3
            diffFlag = true;
        end
        
        Dk = cell(1, length(k_array));
        I  = fluidcell.eye(obj.N, obj.Ntypes);
        
        for i = 1:length(k_array)
            Dk{i} = 1i*k_array(i)*I.*obj.v_BG + diffFlag*k_array(i)^2 * obj.DKernel/2;
        end
        
    end
    
    
    function [dtheta_t, dthetak_t, k_array] = propagateTheta(obj, dtheta, t_array, diffFlag)
        % =================================================================
        % Purpose : Propagates a perturbation on top of a homogeneous
        %           background state. Overloads standard propagation method
        %           in order to include diffusion in equation.
        % Input :   dtheta    -- Initial filling function of perturbation.
        %           t_array   -- Array of evalutation times.
        %           diffFlag  -- if 0, disable diffusion.
        % Output:   dtheta_t  -- Matrix containing perturbation at times
        %                           specified in t_array.
        %           dthetak_t -- Matrix containing momentum modes of
        %                           perturbation at times in t_array.
        %           k_array   -- Array containing the mode numbers.
        % =================================================================
                
        if nargin < 4
            diffFlag = true;
        end
        
        
        % Fourier transform the perturbation to get mode composition.
        phi     = fft( double(dtheta), [], 2 );
        phi     = fftshift(phi,2);
        q_array = 1/(obj.x_grid(2)-obj.x_grid(1)) *(1:(obj.M))/obj.M;
        q_array = q_array - mean(q_array);
        
        k_array = 2*pi*q_array; % array of modes
       
        % Calculate diffusion operator for each mode     
        DK      = obj.calcDiffusiveVelocity( k_array, diffFlag);
        
        f       = zeros(obj.N, obj.N, length(k_array) );
        z       = zeros(1, obj.N, length(k_array) );
        c       = zeros(1, obj.N, length(k_array) );
        
        for k = 1:length(k_array)
            % Get eigenvalues and vectors of diff. operator
            DK_sq       = permute( double(DK{k}), [1 4 2 3 5] );
            [vecs, vals]= eig(DK_sq);
            
            f(:,:,k)    = vecs;
            z(:,:,k)    = diag(vals);
            c(:,:,k)    = vecs\phi(:,k); 
        end
              
        
        % Propagate each mode independently 
        dtheta_t  = zeros(obj.N, obj.M, length(t_array));
        dthetak_t = zeros( length(t_array), obj.M );
        
        for i = 1:length(t_array)            
            phi_t           = sum( f.*c.*exp(-z*t_array(i)) , 2 );
            phi_t           = permute( phi_t , [1 3 2] );
            phi_t           = ifftshift(phi_t, 2);

            dtheta_t(:,:,i) = ifft( phi_t , [], 2 );
            dthetak_t(i,:)  = sum(phi_t.*obj.rapid_w,1);
        end
    end
    

    
end % end public methods

end % end classdef