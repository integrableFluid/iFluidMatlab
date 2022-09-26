classdef Quasi1D_CollisionIntegral < AdvectionSource
    
    
properties (Access = protected)

    lperp   = []; % transverse oscillator length
    gamma   = []; % heating rate
    nu_init = []; % initial transverse populations
    
    gridmap_p = [];
    gridmap_m = [];
    Pp = [];
    Pm = [];
    
end % end protected properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = Quasi1D_CollisionIntegral(model, lperp, varargin)        
        obj = obj@AdvectionSource(model, varargin{:});
             
        obj.lperp = lperp;
        
        % Parse inputs from varargin
        parser = inputParser;
        parser.KeepUnmatched = true;

        addParameter(parser, 'gamma', 0, @(x) x >= 0);
        addParameter(parser, 'nu_init', [0, 0], @(x) all(x >= 0));
        parse(parser,varargin{:});
        
        % add parser results to settings
        obj.gamma   = parser.Results.gamma;
        obj.nu_init = parser.Results.nu_init;
        
        
       % calculate collision grids
        rapid_aux = permute(obj.rapid_grid, [4 2 3 1]);
        k       = 0.5*abs(obj.rapid_grid - rapid_aux);
        qp      = sqrt(k.^2 + 2*lperp^(-2) );
        qm      = real( sqrt(k.^2 - 2*lperp^(-2) ) );
        
        rapid_p = 0.5*(obj.rapid_grid + rapid_aux) + sign(obj.rapid_grid - rapid_aux).*qp;
        rapid_m = 0.5*(obj.rapid_grid + rapid_aux) + sign(obj.rapid_grid - rapid_aux).*qm;
    
        % setup gridmap for interpolation
        obj.gridmap_p = obj.calcInterpolationMap( obj.rapid_grid , rapid_p(:) , 0 , 0);
        obj.gridmap_m = obj.calcInterpolationMap( obj.rapid_grid , rapid_m(:) , 0 , 0);
        
        % calculate common factors in integrals 
        % NOTE: ASSUMES CONSTANT COUPLINGS!!!
        Pp      = obj.model.calcExcitationProb(0, obj.x_grid, k, qp);
        Pm      = obj.model.calcExcitationProb(0, obj.x_grid, k, qm);
        
        obj.Pm = 2 * 2*pi * Pm.*k.*heaviside( 2*k*lperp-2*sqrt(2) );
        obj.Pp = 2 * 2*pi * Pp.*k;
    end
    

    function [I, J, Ip_minus, Ih_minus, Ip_plus, Ih_plus] = calcCollisionIntegral(obj, theta, nu, t)
        % =================================================================
        % Purpose : Calculate collision integral in LL model
        % Input :   theta-- filling function
        %           nu   -- excistation probability (scalar)
        %           t    -- time
        % Output:   I    -- advection collision integral
        %           J    -- excitation "collision integral"
        %           I... -- collision integral for individual channels
        % =================================================================

        [rhoP, rhoS] = obj.model.transform2rho(theta, t);
        
        % calculate rho_p and rho_h on collision grids        
        theta_p   = obj.interp2map( theta, obj.gridmap_p);
        theta_m   = obj.interp2map( theta, obj.gridmap_m);        
        
        % Calculate colition integral components
        Ip_minus = trapz( obj.rapid_grid, double(obj.Pm.*rhoS.t().*(theta.*theta.t().*(1-theta_m).*(1-theta_m.t()))), 4);
        Ih_minus = trapz( obj.rapid_grid, double(obj.Pm.*rhoS.t().*((1-theta).*(1-theta.t()).*theta_m.*theta_m.t())), 4);
        Ip_plus  = trapz( obj.rapid_grid, double(obj.Pp.*rhoS.t().*(theta.*theta.t().*(1-theta_p).*(1-theta_p.t()))), 4);
        Ih_plus  = trapz( obj.rapid_grid, double(obj.Pp.*rhoS.t().*((1-theta).*(1-theta.t()).*theta_p.*theta_p.t())), 4);

        % Normalize minus grids to plus grids
        Nh_plus  = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ih_plus.*rhoS), 1), 2);
        Np_plus  = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ip_plus.*rhoS), 1), 2);
        Nh_minus = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ih_minus.*rhoS), 1), 2);
        Np_minus = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ip_minus.*rhoS), 1), 2);

        Ip_minus = (Nh_plus./Np_minus).*Ip_minus;
        Ih_minus = (Np_plus./Nh_minus).*Ih_minus;
        
        
        N_ext    = [2, 1]; % number of excited atoms in a collision
        zeta     = [0.5, 0.5]; % relative transition strength
        I        = zeros(size(rhoP));
        J        = zeros(size(nu));
        for i = 1:2           
            % Add contributions to collision integral
            Nat         = trapz(obj.x_grid, trapz(obj.rapid_grid, double(rhoP), 1), 2);
            I           = I + zeta(i)*(- Ip_minus + Ih_minus.*(nu(i).^N_ext(i)) - Ip_plus.*(nu(i).^N_ext(i)) + Ih_plus);
            J_temp      = N_ext(i)/(2*Nat) * ( zeta(i)*Ih_plus.*rhoS - zeta(i)*Ip_plus.*rhoS.*(nu(i).^N_ext(i)) );
            J(i)        = trapz(obj.x_grid,  trapz(obj.rapid_grid, double(J_temp), 1), 2);
        end
        
    end


    % Implementation of abstract functions
    
    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        % =================================================================
        % Purpose : Calculates and stores all quantities needed for the
        %           step() method.
        %           In this case of first order step, nothing is required.
        % Input :   theta_init -- Initial filling function.
        %           u_init     -- Initial auxiliary variable 1.
        %           w_init     -- Initial auxiliary variable 2.
        %           t_array    -- Array of time steps.
        % Output:   theta      -- Input theta for first step().
        %           u          -- Input u for first step().
        %           w          -- Input w for first step().
        % =================================================================
       
        if iscell(theta_init) && length(theta_init) > 1
            % Assume theta_init{end} = theta(t = 0),
            % while theta_init{end-1} = theta(t = -dt), ...
            dt          = t_array(2) - t_array(1);
            
            [S, A]      = obj.calcSourceTerm(theta_init{end-1}, obj.nu_init, w_init, -dt);
            
            obj.S_prev  = S;
            obj.A_prev  = A;
            
            theta       = theta_init{end};
            u           = obj.nu_init; 
            w           = w_init; 
        else
            
            [S, A]      = obj.calcSourceTerm(theta_init, obj.nu_init, w_init, 0);
            
            obj.S_prev  = S;
            obj.A_prev  = A;
            
            theta       = theta_init;
            u           = obj.nu_init; 
            w           = w_init;
        end

    end
    
    
    function [S, A] = calcSourceTerm(obj, theta, u, w, t)
        % u is excitation probability nu
        [S, A] = obj.calcCollisionIntegral(theta, u, t);
    end
      
    
    function [u_next, w_next] = auxStep(obj, u, w, A, t, dt)
        % Propagates the excitation probability nu (here u)
        
        u_next = u + dt*A + dt*[obj.gamma*(1-u(1)), obj.gamma*u(1)];
        w_next = w; % second auxiliary variable not used
    end
    
       
end 


methods (Access = private)
    

    function Qint = interp2map(obj, Q, map)
        % =================================================================
        % Purpose : Takes a quantity Q defined on rapidity grid and
        %           maps it to another grid via matrix multiplication.
        % Input :   Q       -- Quantity on rapidity grid.
        %           map     -- (N^2 , N) matrix encoding interpolation.
        % Output:   Qint    -- Quantity on new grid.
        % =================================================================
        
        Q       = double(Q);
        M       = size(Q, 2);

        Q_int = zeros(obj.N, obj.N, M);
        
        for i = 1:M
            Q_temp = map*Q(:,i);
            Q_int(:,:,i) = reshape(Q_temp, obj.N, obj.N);
        end
        
        Qint = fluidcell(permute(Q_int, [1 3 4 2]));
    end
    
    
    function IM = calcInterpolationMap(obj, grid_from, grid_to, extrapflag, cubicflag)
        % =================================================================
        % Purpose : Create map between two grids 
        % Input :   grid_from   -- Old grid.
        %           grid_to     -- New grid.
        %           cubicflag   -- if 1, use cubic splines, else linear
        %                           splines + sparse matrix.
        % Output:   IM           -- Matrix encoding the mapping.
        % =================================================================

        N1 = length(grid_from);
        N2 = length(grid_to);

        IM = zeros(N2,N1);   


        h = diff(grid_from);
        h(end+1) = h(end);

        % Compute "Hessian" matrix
        if cubicflag
            hp = circshift(h,-1);
            hp1 = hp.*(hp+h);
            hp2 = h.*(hp+h);

            lambda = hp./(h+hp);
            mu = 1 - lambda;

            A = diag(2*ones(1,N1)) + diag(lambda(1:end-1),1)  + diag(mu(2:end),-1);
            D = 6*diag( -1./(hp.*h) ) + 6*diag(1./( hp1(1:end-1) ), 1) + 6*diag(1./( hp2(2:end) ), -1);

            B = inv(A)*D;

            B(1,:) = 0;
            B(end,:) = 0;
        else
            B = zeros(N1,N1);
        end

        % Find spline for each point
        for i = 1:N2
            S = cubicflag;

            [ ~, ix ] = min( abs( grid_from-grid_to(i) ) );

            if ix == 1 % extrapolation
                if extrapflag
                    ix = ix+1;
                    S = 0;
                else
                   IM(i,:) = 0;
                   continue
                end
            elseif ix == N1
                if extrapflag
                    ix = ix-1;
                    S = 0;
                else
                   IM(i,:) = 0;
                   continue
                end
            end

            if grid_from(ix) < grid_to(i) % align between gridpoints x_{i-1} and x_{i}
                ix = ix+1;
            end


            dx1 = grid_from(ix) - grid_to(i);
            dx2 = grid_to(i) - grid_from(ix-1);

            IM(i,ix-1) = dx1/h(ix);
            IM(i,ix) = 1-dx1/h(ix);


            IM(i,:) = IM(i,:) + S*B(ix-1,:)*dx1^3 /(6*h(ix));
            IM(i,:) = IM(i,:) + S*B(ix,:)*dx2^3 /(6*h(ix));
            IM(i,:) = IM(i,:) - S*B(ix-1,:)*dx1*h(ix)/6;
            IM(i,:) = IM(i,:) - S*B(ix,:)*dx2*h(ix)/6;
        end

        if ~cubicflag
            IM = sparse(IM);
        end

    end
    
end % end private methods


end % end classdef