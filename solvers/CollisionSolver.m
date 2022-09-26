classdef CollisionSolver < iFluidSolver
    
% Solves GHD equation with collision integral using a split-step
% propagation scheme.
% Currently only works for the Lieb-Liniger model.
    
properties (Access = protected)
    theta_mid = []; % Midpoint filling required for taking 2nd order step
    nu_mid = [];
    
    I_prev = []; 
    J_prev = [];
    I_mid = []; 
    J_mid = [];

    
    lperp = [];
    gamma = [];
    nu_init = [0, 0];
    
    gridmap_p = [];
    gridmap_m = [];
    Pp = [];
    Pm = [];
    
end % end protected properties


methods (Access = public)
    
    % Superclass constructor
    function obj = CollisionSolver(coreObj, lperp, gamma, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
        
        obj.lperp = lperp;
        obj.gamma = gamma;
        obj.calcCharac = false;
                
        
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
        Pp      = obj.coreObj.calcExcitationProb(0, obj.x_grid, k, qp);
        Pm      = obj.coreObj.calcExcitationProb(0, obj.x_grid, k, qm);
        
        obj.Pm = 2*(2*pi)^2 * Pm.*k.*heaviside( 2*k*lperp-2*sqrt(2) );
        obj.Pp = 2*(2*pi)^2 * Pp.*k;
        
    end
   
    
    function setInitialNu(obj, nu_init)
        obj.nu_init = nu_init;        
    end
    
    
    function [I, J, Ip_minus, Ih_minus, Ip_plus, Ih_plus] = calcCollisionIntegral(obj, rhoP, rhoH, nu)
        % =================================================================
        % Purpose : Calculate collision integral in LL model
        % Input :   rhoP -- root density (iFluidTensor)
        %           rhoH -- density of holes (iFluidTensor)
        %           nu   -- excistation probability (scalar)
        % Output:   I    -- quasi-particle collision integral
        %           J    -- excitation "collision integral"
        % =================================================================

        % calculate rho_p and rho_h on collision grids        
        rhoP_m   = obj.interp2map( rhoP, obj.gridmap_m);
        rhoP_p   = obj.interp2map( rhoP, obj.gridmap_p);
        rhoH_p   = obj.interp2map( rhoH, obj.gridmap_p);
        rhoH_m   = obj.interp2map( rhoH, obj.gridmap_m);
        
        % Calculate colition integral components
        Ip_minus = trapz( obj.rapid_grid, double(obj.Pm.*(rhoP.*rhoP.t().*rhoH_m.*rhoH_m.t())), 4);
        Ih_minus = trapz( obj.rapid_grid, double(obj.Pm.*(rhoH.*rhoH.t().*rhoP_m.*rhoP_m.t())), 4);
        Ip_plus  = trapz( obj.rapid_grid, double(obj.Pp.*(rhoP.*rhoP.t().*rhoH_p.*rhoH_p.t())), 4);
        Ih_plus  = trapz( obj.rapid_grid, double(obj.Pp.*(rhoH.*rhoH.t().*rhoP_p.*rhoP_p.t())), 4);

        % Normalize minus grids to plus grids
        Nh_plus  = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ih_plus), 1), 2);
        Np_plus  = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ip_plus), 1), 2);
        Nh_minus = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ih_minus), 1), 2);
        Np_minus = trapz(obj.x_grid, trapz(obj.rapid_grid, double(Ip_minus), 1), 2);

        Ip_minus = (Nh_plus./Np_minus).*Ip_minus;
        Ih_minus = (Np_plus./Nh_minus).*Ih_minus;
        
        
        Next     = [2, 1]; % number of excited atoms in a collision
        zeta     = [0.5, 0.5]; % relative transition strength
        I        = zeros(size(rhoP));
        J        = zeros(size(nu));
        for i = 1:2           
            % Add contributions to collision integral
            Nat         = trapz(obj.x_grid, trapz(obj.rapid_grid, double(rhoP), 1), 2);
            I           = I + zeta(i)*(- Ip_minus + Ih_minus.*(nu(i).^Next(i)) - Ip_plus.*(nu(i).^Next(i)) + Ih_plus);
            J_temp      = Next(i)/(2*Nat) * ( zeta(i)*Ih_plus - zeta(i)*Ip_plus.*(nu(i).^Next(i)) );
            J(i)        = trapz(obj.x_grid,  trapz(obj.rapid_grid, double(J_temp), 1), 2);
        end
        
    end
    
end % end public methods


methods (Access = protected)

    function [theta, u, w] = initialize(obj, theta_init, u_init, w_init, t_array)
        % =================================================================
        % Purpose : Calculates and stores all quantities needed for the
        %           step() method.
        %           In this the filling at dt/2 is needed to start the
        %            propagation.
        % Input :   theta_init -- Initial filling function.
        %           u_init     -- Initial position characteristic.
        %           w_init     -- Initial rapidity characteristic.
        %           t_array    -- Array of time steps.
        % Output:   theta      -- Input theta for first step().
        %           u          -- Input u for first step().
        %           w          -- Input w for first step().
        % =================================================================
        dt      = t_array(2) - t_array(1);
        ddt     = dt/2/10;
        theta   = theta_init;
        u       = obj.nu_init;          % u is excitation probability in this solver
        w       = w_init;

        % Set up collision integrals for propagation
        [rhoP, rhoS]= obj.coreObj.transform2rho(theta_init, 0);
        rhoH        = rhoS - rhoP;
        [I, J]      = obj.calcCollisionIntegral( rhoP, rhoH, u);
        obj.I_prev  = I;
        obj.J_prev  = J;
        obj.I_mid   = I;
        obj.J_mid   = J;
        
        % Calculate first theta_mid at t = dt/10/2 using first order
        % step, then use that to calculate the actual theta_mid at
        % t = dt/2 using second order steps. 
        obj.theta_mid = theta_init + ddt/2*I;
        obj.nu_mid    = u + ddt/2*J + ddt/2*[obj.gamma, obj.gamma*u(1)];      
        obj.theta_mid = obj.performFirstOrderStep(obj.theta_mid, u_init, w_init, 0, ddt/2);
        
        theta_temp  = theta;
        u_temp      = u;
        
        for i = 1:10
            t           = (i-1)*ddt;
            [theta_temp, u_temp] = obj.step(theta_temp, u_temp, w_init, t, ddt);
        end

        obj.theta_mid = theta_temp;
        obj.nu_mid = u_temp;
        
        obj.I_mid = obj.I_prev;
        obj.J_mid = obj.J_prev;
        
        obj.I_prev = I;
        obj.J_prev = J;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_next, t, dt)
        % =================================================================
        % Purpose : Performs a split step propagating
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
        
        % Calculate collision integral
        [rhoP, rhoS]= obj.coreObj.transform2rho(theta_prev, t);
        rhoH        = rhoS - rhoP;
        [I, J]      = obj.calcCollisionIntegral( rhoP, rhoH, u_prev);

        % First part of split step: solve collision equation
        rhoP        = rhoP + 0.5*dt*(3*I - obj.I_prev);
        u_next      = u_prev + 0.5*dt*(3*J - obj.J_prev) + dt*[obj.gamma*(1-u_prev(1)), obj.gamma*u_prev(1)];
        
        obj.I_prev  = I;
        obj.J_prev  = J;
        
        theta_prev  = obj.coreObj.transform2theta(rhoP, t);       
        
        % Second part of split step: solve propagation equation
        theta_next = step2(obj, obj.theta_mid, theta_prev, t, dt);


        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % I.e. calculate theta(t+dt+dt/2) using theta(t+dt) as midpoint.
        [rhoP, rhoS]= obj.coreObj.transform2rho(obj.theta_mid, t+dt/2);
        rhoH        = rhoS - rhoP;
        [I, J]      = obj.calcCollisionIntegral( rhoP, rhoH, obj.nu_mid);

        % First part of split step: solve collision equation
        rhoP        = rhoP + 0.5*dt*(3*I - obj.I_mid);
        obj.nu_mid  = obj.nu_mid + 0.5*dt*(3*J - obj.J_mid) + dt*[obj.gamma*(1-obj.nu_mid(1)), obj.gamma*obj.nu_mid(1)];
        
        obj.I_mid   = I;
        obj.J_mid   = J;
        
        obj.theta_mid  = obj.coreObj.transform2theta(rhoP, t+dt/2);         
        obj.theta_mid  = step2(obj, theta_next, obj.theta_mid, t+dt/2, dt);


        function theta_next = step2(obj, theta_mid, theta_prev, t, dt)
            % Estimate x' and rapid' using midpoint filling
            [v_eff, a_eff] = obj.coreObj.calcEffectiveVelocities(theta_mid, t+dt/2, obj.x_grid, obj.rapid_grid, obj.type_grid);

            x_mid   = obj.x_grid - 0.5*dt*v_eff;
            r_mid   = obj.rapid_grid - 0.5*dt*a_eff; 

            % Interpolate v_eff and a_eff to midpoint coordinates x' and rapid' 
            v_mid   = obj.interpPhaseSpace( v_eff, r_mid, x_mid, true ); % always extrapolate v_eff
            a_mid   = obj.interpPhaseSpace( a_eff, r_mid, x_mid, true );
            
            % From midpoint velocity calculate fullstep x and rapid translation
            x_back  = obj.x_grid - dt*v_mid;
            r_back  = obj.rapid_grid - dt*a_mid;

            % Use interpolation to find theta_prev at x_back, r_back and
            % assign values to theta_next.
            theta_next = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);

        end % end nested function
    end
    
    
end % end protected methods


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