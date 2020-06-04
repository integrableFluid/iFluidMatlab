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
        % Not used
    end
      
    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_prev, t, dt)
        % Not used
    end
    
    
end % end protected methods


methods (Access = public)
    
    function [D, v_eff, wt] = calcDiffusionKernel(obj, theta_BG)
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        [rho_BG, rhoS_BG] = obj.coreObj.transform2rho(theta_BG, 0);
        
        I           = iFluidTensor(obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes, 'eye')./obj.rapid_w;
        T           = -1/(2*pi)*obj.coreObj.getScatteringRapidDeriv(0, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);
%         T = T.*sqrt(obj.rapid_w.*permute(obj.rapid_w, [4 2 3 1]));
        T_dr        = obj.coreObj.applyDressing(T, theta_BG, 0);
        v_eff       = obj.coreObj.calcEffectiveVelocities(theta_BG, 0, obj.x_grid, obj.rapid_grid, obj.type_grid);
        f           = obj.coreObj.getStatFactor( theta_BG );
        
        W           = rho_BG.*f.*T_dr.^2 .* abs(v_eff - v_eff.t());
%         W = W.*sqrt(obj.rapid_w.*permute(obj.rapid_w, [4 2 3 1]));
%         W = W.*obj.rapid_w;
        
        w           = sum( transpose(W.*obj.rapid_w) , 4);
%         w           = sum( W.t() , 4);
        
        D           = rhoS_BG.^(-2).*(I.*w - W).*sqrt(obj.rapid_w.*permute(obj.rapid_w, [4 2 3 1]));
        

        wt          = 0.5*rhoS_BG.^(-2).*w;
%         
%         D           = D.getX(1);
%         v_eff       = v_eff.getX(1);
    end
    
    
    function Dk = calcDiffusiveVelocity(obj, k_array, diffFlag)
        if nargin < 3
            diffFlag = true;
        end
        
        Dk = cell(1, length(k_array));
        
        I       = iFluidTensor(obj.N, 1, obj.Ntypes, obj.N, obj.Ntypes, 'eye');
        
        for i = 1:length(k_array)
            Dk{i} = 1i*k_array(i)*I.*obj.v_BG + diffFlag*k_array(i)^2 * obj.DKernel/2;
        end
        
    end
    
    
    function [dtheta_t, dthetak_t, k_array] = propagateTheta(obj, dtheta, t_array, diffFlag)
        % OVERLOAD
        
        %  *** First find coefficients c ***
        % Fourier transform dtheta to get g
        % Calc eigenfunctions of Dk for each relevant k --> f
        % Invert matrix of f on k and rapid (or better use \) to find c's
        
        % *** perform time evolution ***
        
        if nargin < 4
            diffFlag = true;
        end
        
        
        %% Fourier stuff here
        phi = fft( double(dtheta), [], 2 );
        P1 = phi;
        P1 = fftshift(phi,2);
        q_array = 1/(obj.x_grid(2)-obj.x_grid(1)) *(1:(obj.M))/obj.M;
        q_array = q_array - mean(q_array);
        
        k_array = 2*pi*q_array;
        
        figure
        subplot(1,2,1)
        plot(k_array, real(P1(obj.N/2,:)))
        subplot(1,2,2)
        plot(k_array, imag(P1(obj.N/2,:)))

        
        %% Eigen stuff here        
        DK = obj.calcDiffusiveVelocity( k_array, diffFlag);
        
        f = zeros(obj.N, obj.N, length(k_array) );
        z = zeros(1, obj.N, length(k_array) );
        c = zeros(1, obj.N, length(k_array) );
        
        d = zeros(1, obj.N, length(k_array) );
        
        for k = 1:length(k_array)
            DK_sq = permute( double(DK{k}), [1 4 2 3 5] );
            [vecs, vals] = eig(DK_sq);
            
            % find orthonormal eigenbasis and transformaiton matrix
%             g = orth(vecs);
%             At = g\vecs;
%             A = transpose(At); % f^T = A*g^T
%             ATi = inv(At);
%             
%             % find b coefficients --> d coefficients
% %             b = (vecs*ATi)\P1(:,k);
%             b = g\P1(:,k);
%             d(:,:,k) = ATi*b;
% %             d(:,:,k) = transpose(b)*inv(A);


%             g = orth(vecs);
%             AT = g' * vecs;
%             b = g' * P1(:,k);
%             d(:,:,k) = AT\b;

            
            f(:,:,k) = vecs;
            z(:,:,k) = diag(vals);
            c(:,:,k) = vecs\P1(:,k); 
        end
        
        

        %% Time evolution here
%         dtheta_t = zeros(obj.N, obj.M, length(t_array));
%         dthetak_t = zeros( length(t_array), obj.M );
%         for i = 1:length(t_array)
%                         
%             phi_t       = sum( f.*d.*exp(-z*t_array(i)) , 2 );
%             phi_t       = permute( phi_t , [1 3 2] );
%             
% %             phi_t       = phi_t.*obj.rapid_w.^(-1/2);
%             phi_t = ifftshift(phi_t, 2);
%             
%             dtheta_t(:,:,i) = ifft( phi_t , [], 2 );
%             dthetak_t(i,:)= sum(phi_t.*obj.rapid_w,1);
%         end
        
%         figure
%         imagesc( squeeze(real(d-c)) )
%         figure
%         imagesc( squeeze(imag(d-c)) )
        
        
        
        
        dtheta_t = zeros(obj.N, obj.M, length(t_array));
        dthetak_t = zeros( length(t_array), obj.M );
        
        for i = 1:length(t_array)
                        
            phi_t       = sum( f.*c.*exp(-z*t_array(i)) , 2 );
            phi_t       = permute( phi_t , [1 3 2] );
            
            phi_t = ifftshift(phi_t, 2);

            dtheta_t(:,:,i) = ifft( phi_t , [], 2 );
            dthetak_t(i,:)= sum(phi_t.*obj.rapid_w,1);
        end
    end
    

    
end % end public methods

end % end classdef