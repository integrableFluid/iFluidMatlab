classdef CollisionSolverTEST < iFluidSolver
    
properties (Access = protected)
    theta_mid = []; % Midpoint filling required for taking 2nd order step
    lperp = [];
    gridmap_p = [];
    gridmap_m = [];
    P1 = [];
    P2 = [];
    
    P3 = [];
    P4 = [];
    
end % end protected properties


methods (Access = public)
    
    % Superclass constructor
    function obj = CollisionSolverTEST(coreObj, lperp, Options)        
        obj = obj@iFluidSolver(coreObj, Options);
        
        obj.lperp = lperp;
        obj.calcCharac = false;
                
        
        % calculate collision grids
        rapid_aux = permute(obj.rapid_grid, [4 2 3 1]);
        k       = 0.5*abs(obj.rapid_grid - rapid_aux);
        qp      = sqrt(k.^2 + 2*lperp^(-2) );
        qm      = real( sqrt(k.^2 - 2*lperp^(-2) ) );
        
        rapid_p = 0.5*(obj.rapid_grid + rapid_aux) + sign(obj.rapid_grid - rapid_aux).*qp;
        rapid_m = 0.5*(obj.rapid_grid + rapid_aux) + sign(obj.rapid_grid - rapid_aux).*qm;
    
        % setup gridmap for interpolation
        obj.gridmap_p = obj.calcInterpolationMap( obj.rapid_grid , rapid_p(:) , 1);
        obj.gridmap_m = obj.calcInterpolationMap( obj.rapid_grid , rapid_m(:) , 1);
        
        % calculate common factors in integrals 
        % NOTE: ASSUMES CONSTANT COUPLINGS!!!
        Pp      = obj.coreObj.calcExctProb(0, obj.x_grid, k, qp);
        Pm      = obj.coreObj.calcExctProb(0, obj.x_grid, k, qm);
        
        obj.P1 = 2*(2*pi)^2 * Pm.*2.*k.*heaviside( 2*k*lperp-2*sqrt(2) );
        obj.P2 = 2*(2*pi)^2 * Pp.*2.*k;
        
        
        
        obj.P3 = (2*pi)^2 * Pp .* (2*k);
        couplings = obj.coreObj.getCouplings();
        c = couplings{1,2}(0,0);
        c = obj.coreObj.convert2SI(c,'length');
        obj.P4 = 4*pi*c^2 *2*k./( c^2 + (obj.rapid_grid - rapid_aux).^2 );

    end
    
    
    function [I, J] = calcCollisionIntegral(obj, t, theta, nu)

        % calculate rho_p and rho_h on collision grids
        [rhoP, rhoS] = obj.coreObj.transform2rho(theta, t);
        rhoH    = rhoS - rhoP;
        
        rhoP_m  = obj.interp2map( rhoP, obj.gridmap_m);
        rhoP_p  = obj.interp2map( rhoP, obj.gridmap_p);
        rhoH_p  = obj.interp2map( rhoH, obj.gridmap_p);
        rhoH_m  = obj.interp2map( rhoH, obj.gridmap_m);
        
        nu_m    = obj.interp2map( nu, obj.gridmap_m);

        % evaluate collision integral
        Iint    = obj.P1.*(-theta.*rhoP.t().*rhoH_m.*rhoH_m.t()        + ...
                     0.5*(1-theta).*rhoH.t().*rhoP_m.*rhoP_m.t().*(nu_m+nu_m.t()))+ ...
                  obj.P2.*(-0.5*theta.*rhoP.t().*rhoH_p.*rhoH_p.t().*(nu+nu.t())    + ...
                            (1-theta).*rhoH.t().*rhoP_p.*rhoP_p.t()         );
                   
        % integrate Iint over rapid_aux
        I       = sum(Iint.*permute(obj.rapid_w, [4 2 3 1]), 4);
        
        % calculate collision integral for nu   
        Jint    = obj.P3.*rhoP_p.*rhoP_p.t().*rhoH.*rhoH.t() - ...
                  obj.P3.*rhoH_p.*rhoH_p.t().*rhoP.*rhoP.t().*nu + ...
                  obj.P4.*rhoP.t().*( nu.t() - nu );                
                                   
        J       = sum(Jint.*permute(obj.rapid_w, [4 2 3 1]), 4);
        
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
        u       = zeros( size(theta) );
        u       = iFluidTensor(u);
        w       = w_init;

        % Calculate first theta_mid at t = dt/10/2 using first order
        % step, then use that to calculate the actual theta_mid at
        % t = dt/2 using second order steps. 
        obj.theta_mid = obj.performFirstOrderStep(theta_init, u_init, w_init, 0, ddt);
        theta_temp = theta;
        u_temp = u;

        for i = 1:10
            t           = (i-1)*ddt;
            [theta_temp, u_temp] = obj.step(theta_temp, u_temp, w_init, t, ddt);
        end

        obj.theta_mid = theta_temp;
    end
      

    function [theta_next, u_next, w_next] = step(obj, theta_prev, u_prev, w_next, t, dt)
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
        
        
        % Split step
        [I, J] = obj.calcCollisionIntegral( t, theta_prev, u_prev);
        
        u_prev = J*dt + u_prev;
        theta_prev = I*dt + theta_prev;
        
        
        % First, calculate theta(t+dt) using the midpoint filling theta(t+dt/2)
        [theta_next, u_next] = step2(obj, obj.theta_mid, theta_prev, u_prev, t, dt);
        
        
        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % I.e. calculate theta(t+dt+dt/2) using theta(t+dt) as midpoint.
        obj.theta_mid  = step2(obj, theta_next, obj.theta_mid, u_next, t+dt/2, dt);
        
        
        function [theta_next, nu_next] = step2(obj, theta_mid, theta_prev, nu_prev, t, dt)
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
            nu_next = obj.interpPhaseSpace(nu_prev, r_back, x_back, true);
        end % end nested function
    end
    
    
end % end protected methods


methods (Access = private)
    
    function Qint = interp2map(obj, Q, map)
        % Q exist on rapid_grid
        Q = double(Q);
        M = size(Q, 2);
        
        Q_int = zeros(obj.N, obj.N, M);
        
        for i = 1:M
            Q_temp = map*Q(:,i);
            Q_int(:,:,i) = reshape(Q_temp, obj.N, obj.N);
        end
        
        Qint = iFluidTensor(permute(Q_int, [1 3 4 2]));
    end
    
    function IM = calcInterpolationMap( obj, grid, array , cubicflag)
        N1 = length(grid);
        N2 = length(array);

        IM = zeros(N2,N1);   


        h = diff(grid);
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

            [ ~, ix ] = min( abs( grid-array(i) ) );

            if ix == 1 % extrapolation
                ix = ix+1;
                S = 0;
            elseif ix == N1
                ix = ix-1;
                S = 0;
            end

            if grid(ix) < array(i) % align between gridpoints x_{i-1} and x_{i}
                ix = ix+1;
            end


            dx1 = grid(ix) - array(i);
            dx2 = array(i) - grid(ix-1);

            IM(i,ix-1) = dx1/h(ix);
            IM(i,ix) = 1-dx1/h(ix);


            IM(i,:) = IM(i,:) + S*B(ix-1,:)*dx1^3 /(6*h(ix));
            IM(i,:) = IM(i,:) + S*B(ix,:)*dx2^3 /(6*h(ix));
            IM(i,:) = IM(i,:) - S*B(ix-1,:)*dx1*h(ix)/6;
            IM(i,:) = IM(i,:) - S*B(ix,:)*dx2*h(ix)/6;
        end

    end
    
end % end private methods

end % end classdef