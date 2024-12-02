classdef SineGordonModel < iFluidCore
    % iFluid implmentation of TBA functions for sine-Gordon model
    %
    % ### Couplings are {mu}
    %
    % Quasiparticle species are indexed in the following order:
    %   [Breathers, Soliton, level-1 magnons, level-2 magnons, ...]
    %

properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies    = 'fermion'; 
    
    N_bre           = [];  % number of breather types
    N_mag           = [];  % number of magnon types
    M_sol           = 1;   % mass of soliton 

    alpha           = []; % inverse continuous fraction of breathers 
    xi              = []; % renormalized coupling

    l_var           = [];
    v_var           = [];
    r_var           = [];

    scat_kern       = []; % scattering kernel
    per_scat_kern   = []; % periodic scattering kernel
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = SineGordonModel(x_grid, rapid_grid, rapid_w, particles, couplings, Options)   
        % Set default parameters
        if nargin < 5
            Options = struct;
        end
        
        % The array particles specifies the number of each quasi-particle
        % type. Should be structures as:
        %   particles = [N_bre, N_mag1, N_mag2, ... ]
        assert(length(particles) >= 2, 'Must specify at least number of breathers and level-1 magnons!')

        N_bre   = particles(1);
        N_mag   = particles(2:end);
        Ntypes  = sum(N_mag) + N_bre + 1; % plus 1 soliton type 

        % Call superclass constructor
        obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
        
        % Calculate and store coupling
        alpha   = Inf;
        for i = length(N_mag):-1:1
            alpha   = N_mag(i) + 1/alpha;
        end

        obj.N_mag = N_mag;
        obj.N_bre = N_bre;

        obj.alpha = alpha;
        obj.xi    = 1/( N_bre + 1/alpha );

        [obj.l_var, obj.v_var, obj.r_var] = obj.getAuxVariablesTBA();
       
        obj.scat_kern = fluidcell( obj.getTotalScatteringKernel() );
        obj.per_scat_kern = fluidcell( obj.getTotalPeriodicScatteringKernel() );
    end
    
    
    %% Implementation of abstract equations
    
    function M = getMasses(obj, type)
        % get mass of breather/soliton with the given type index
        % type = 1:N_bre --> breathers
        % type = N_bre + 1 --> soliton
        % type > N_bre + 1 -- > magnons
  
        M = 2*obj.M_sol*sin(type*pi*obj.xi/2); % breather masses
        
        M(type == obj.N_bre + 1) = obj.M_sol; % soliton mass
        M(type > obj.N_bre + 1) = 0; % magnon masses
    end


    function xi = getXi(obj)
        xi = obj.xi;
    end


    function [l, v, r, y, p, kappa] = getAuxVariablesTBA(obj)
        % Returns additional variables needed for the fully coupled TBA 
        % equations.
 
        level = length(obj.N_mag);
        
        p = [obj.alpha, 1];
        for i = 1:level
            p = [p, p(end-1) - p(end) * floor(p(end-1)/p(end))];
        end
        
        y = [0, 1];
        if level > 0
            y = [y, obj.N_mag(1)];
        end
        for i = 4:level+2
            y = [y, y(end-1) + obj.N_mag(i-2) * y(end)];
        end
        
        kappa = [0, cumsum(obj.N_mag)];
        
        l = zeros(1, sum(obj.N_mag));
        v = zeros(1, sum(obj.N_mag));
        ii = 0;
        for jj = 1:sum(obj.N_mag)
            l(jj) = y(ii + 2 - 1) + (jj - kappa(ii + 1)) * y(ii + 2);
            v(jj) = (-1)^(floor((l(jj)-1)/p(1)));
            if jj == kappa(ii + 1 + 1)
                l(jj) = y(ii + 2 + 1 -1);
                v(jj) = (-1)^(ii+1);
                ii = ii + 1;
            end
        end
        
        r = zeros(1, sum(obj.N_mag));
        ii = 0;
        for jj = 1:sum(obj.N_mag)
            if jj == sum(obj.N_mag) || jj == kappa(ii + 1 + 1)
                ii = ii + 1;
            end
            r(jj) = (-1)^ii * (p(ii + 1) - (jj - kappa(ii + 1)) * p(ii + 1 + 1));
        end



    end


    function eta = getTypeSigns(obj)
        % Returns parity of each quasiparticle species
        eta                 = zeros(1,1,obj.Ntypes);
        eta(1:obj.N_bre)    = 1; % breathers
        eta(obj.N_bre + 1)  = 1; % soliton
        eta(obj.N_bre+2:end)= -sign(obj.r_var); % magnons
        
        eta                 = fluidcell(eta);
    end


    function qbare = getBareTopologicalCharge(obj)
        qbare               = zeros(1,1,obj.Ntypes); % breathers have no charge
        qbare(obj.N_bre+1)  = 1; % soliton
        qbare(obj.N_bre+2:end) = -2*obj.l_var; % magnons
    end
    

    function ebare = getBareEnergy(obj, t, x, rapid, type)
        masses  = obj.getMasses(type);
        mu      = obj.couplings{1,1}(t,x);
        q       = obj.getBareTopologicalCharge();

        ebare   = masses.*cosh(rapid) - q*mu;
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        masses  = obj.getMasses(type);        
        pbare   = masses.*sinh(rapid);
    end
    
    
    function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
        masses  = obj.getMasses(type);        
        de      = masses.*sinh(rapid);
    end

    
    function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
        masses  = obj.getMasses(type);        
        dp      = masses.*cosh(rapid);
    end


    function Phi0 = getScatteringKernel_SS(obj)
        % Returns Soliton-Soliton scattering kernel as N x N matrix
        
        rmax    = -obj.rapid_grid(1);
        
        dr_grid = linspace(-2*rmax, 2*rmax, 2*obj.N+1); % vector containing all possible rapid differences with given grids
        
        T       = pi/(dr_grid(2)-dr_grid(1)); % assuming linear rapidity grid
        t_grid  = linspace(-T, T, 2*obj.N+1); % Fourier grid
        t_grid(end) = [];
        
        % calculate self-coupling as function of rapidity 
        Phi_self_ft = -sinh((1-obj.xi)*pi*t_grid/2)./( 2*sinh(pi*obj.xi*t_grid/2).*cosh(pi*t_grid/2) + eps );
        Phi_self = real( 2*T*fftshift(ifft( Phi_self_ft )).*exp(-1i*dr_grid(1:end-1)*T) );
        
        Phi_self = Phi_self - Phi_self(end); % NOTE: somehow the FT'ed kernel is offset from 0, this fixes it
        % Phi_self(end+1) = 0;

        % map to convolution kernel (rapid_grid - rapid_grid')
        Phi0 = diag( Phi_self(obj.N+1)*ones(1,obj.N)); 
        for i = 2:obj.N
            diag_elem = diag( Phi_self(obj.N+i)*ones(1,obj.N-i+1), i-1) + diag( Phi_self(obj.N-i+2)*ones(1,obj.N-i+1), -i+1);
            Phi0 = Phi0 + diag_elem;
        end

       
    end


    function Phi = getScatteringKernel_SB(obj, i)
        % Returns Soliton-Breather scattering kernel as N x N matrix
        % i = 1...N_bre is breather index

        rapid   = obj.rapid_grid - obj.rapid_grid';
        
        Phi     = -4*cos(i*pi*obj.xi/2)*cosh(rapid)./(cos(i*pi*obj.xi) + cosh(2*rapid));

        for k = 1:i-1
            Phi = Phi - 2*cos((i-2*k)*pi*obj.xi/2)./(cosh(rapid) - sin((i-2*k)*pi*obj.xi/2));
        end
    end


    function Phi = getScatteringKernel_BB(obj, i, j)
        % Returns Breather-Breather scattering kernel as N x N matrix
        % i,j = 1...N_bre are breather index

        rapid = obj.rapid_grid - obj.rapid_grid';

        Phi = 4*sin((i+j)*pi*obj.xi/2)*cosh(rapid)./(cos((i+j)*pi*obj.xi) - cosh(2*rapid) + eps) + ...
              4*sin((i-j)*pi*obj.xi/2)*cosh(rapid)./(cos((i-j)*pi*obj.xi) - cosh(2*rapid) + eps);

        for k = 1:j-1
            Phi = Phi - 2*sin((j-i-2*k)*pi*obj.xi/2)./(cos((j-i-2*k)*pi*obj.xi/2) - cosh(rapid) + eps) - ...
                        2*sin((j+i-2*k)*pi*obj.xi/2)./(cos((j+i-2*k)*pi*obj.xi/2) + cosh(rapid) + eps);
        end
    end


    function a = getScatteringA(obj, rapid, k, p)
        % a-function
        assert(abs(p) == 1) % p is a parity

        if obj.alpha > 0 && round(k/obj.alpha) == k/obj.alpha % is integer
            a = zeros(size(rapid));
        else
            a = pi/obj.alpha*sin(pi*k/obj.alpha)./(cos(pi*k/obj.alpha) - p*cosh(pi*rapid/obj.alpha));
        end
    end


    function Phi = getScatteringKernel_SM(obj, i)
        % Returns Soliton-Magnon scattering kernel as N x N matrix
        % i is magnon index

        rapid   = obj.rapid_grid - obj.rapid_grid';
        chi     = 2*obj.alpha/pi/obj.xi; 

        Phi     = chi*obj.getScatteringA(chi*rapid, obj.l_var(i), obj.v_var(i));
    end


    function Phi = getScatteringKernel_MM(obj, i, j)
        % Returns Magnon-Magnon scattering kernel as N x N matrix
        % i and j are magnon indices

        rapid   = obj.rapid_grid - obj.rapid_grid';
        chi     = 2*obj.alpha/pi/obj.xi; 

        temp    = zeros(obj.N, obj.N);
        kmax    = min([obj.l_var(i), obj.l_var(j)]) - 1;
        for k = 1:kmax
            temp = temp + 2*obj.getScatteringA(chi*rapid, abs(obj.l_var(i)-obj.l_var(j)) + 2*k, obj.v_var(i)*obj.v_var(j));
        end

        Phi     = -chi*obj.getScatteringA(chi*rapid, abs(obj.l_var(i)-obj.l_var(j)), obj.v_var(i)*obj.v_var(j)) + ...
                  -chi*obj.getScatteringA(chi*rapid, obj.l_var(i)+obj.l_var(j), obj.v_var(i)*obj.v_var(j)) + ...
                  -chi*temp;
    end


    function Phi_tot = getTotalScatteringKernel(obj)
        % Returns total scattering kernel

        Phi_tot = zeros(obj.N, obj.N, obj.Ntypes, obj.Ntypes);

        % calculate soliton-soliton scattering (assume for now that this is
        % given by the self-scattering \Phi_0)
        Phi     = obj.getScatteringKernel_SS();
        Phi_tot(:, :, obj.N_bre+1, obj.N_bre+1) = Phi;

        % calculate breather-breather scattering
        for i = 1:obj.N_bre
        for j = 1:obj.N_bre
            Phi     = obj.getScatteringKernel_BB(i, j);
            Phi_tot(:, :, i, j) = Phi;
        end
        end

        % calculate breather-soliton scattering
        for i = 1:obj.N_bre
            Phi     = obj.getScatteringKernel_SB(i);
            Phi_tot(:, :, obj.N_bre+1, i) = Phi;
            Phi_tot(:, :, i, obj.N_bre+1) = Phi;
        end

        % calculate soliton-magnon scattering
        for i = 1:sum(obj.N_mag)
            Phi     = obj.getScatteringKernel_SM(i);
            Phi_tot(:, :, obj.N_bre+1, obj.N_bre+1+i) = Phi;
            Phi_tot(:, :, obj.N_bre+1+i, obj.N_bre+1) = Phi;
        end

        % calculate magnon-magnon scattering
        for i = 1:sum(obj.N_mag)
        for j = 1:sum(obj.N_mag)
            Phi     = obj.getScatteringKernel_MM(i, j);
            Phi_tot(:, :, obj.N_bre+1+i, obj.N_bre+1+j) = Phi;
        end
        end

        Phi_tot = fluidcell(permute(Phi_tot, [1 5 3 2 4]));

    end


    function Phi_tot = getTotalPeriodicScatteringKernel(obj)
        % Returns total scattering kernel, with periodic boundary rapidity
        % boundary conditions for magnon scattering.

        Phi_tot = zeros(obj.N, obj.N, obj.Ntypes, obj.Ntypes);

        % calculate soliton-soliton scattering (assume for now that this is
        % given by the self-scattering \Phi_0)
        Phi     = obj.getScatteringKernel_SS();
        Phi_tot(:, :, obj.N_bre+1, obj.N_bre+1) = Phi;

        % calculate breather-breather scattering
        for i = 1:obj.N_bre
        for j = 1:obj.N_bre
            Phi     = obj.getScatteringKernel_BB(i, j);
            Phi_tot(:, :, i, j) = Phi;
        end
        end

        % calculate breather-soliton scattering
        for i = 1:obj.N_bre
            Phi     = obj.getScatteringKernel_SB(i);
            Phi_tot(:, :, obj.N_bre+1, i) = Phi;
            Phi_tot(:, :, i, obj.N_bre+1) = Phi;
        end

        % calculate soliton-magnon scattering
        for i = 1:sum(obj.N_mag)
            Phi     = obj.getScatteringKernel_SM(i);

            Phi_per = fftshift(Phi);
            Phi_per(1:floor(obj.N/2), 1:floor(obj.N/2)) = 0;
            Phi_per(ceil(obj.N/2):end, ceil(obj.N/2):end) = 0;
            Phi = fluidcell(Phi + Phi_per);


            Phi_tot(:, :, obj.N_bre+1, obj.N_bre+1+i) = Phi;
            Phi_tot(:, :, obj.N_bre+1+i, obj.N_bre+1) = Phi;
        end

        % calculate magnon-magnon scattering
        for i = 1:sum(obj.N_mag)
        for j = 1:sum(obj.N_mag)
            Phi     = obj.getScatteringKernel_MM(i, j);

            Phi_per = fftshift(Phi);
            Phi_per(1:floor(obj.N/2), 1:floor(obj.N/2)) = 0;
            Phi_per(ceil(obj.N/2):end, ceil(obj.N/2):end) = 0;
            Phi = fluidcell(Phi + Phi_per);

            Phi_tot(:, :, obj.N_bre+1+i, obj.N_bre+1+j) = Phi;
        end
        end

        Phi_tot = fluidcell(permute(Phi_tot, [1 5 3 2 4]));

    end
    

    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % For now, scattering kernel is pre-calculated, i.e. rapidity
        % dependence is set by default grid and there is no time or
        % position dependce.
        dT = obj.scat_kern; 
    end


    function dT = getPeriodicScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % For now, scattering kernel is pre-calculated, i.e. rapidity
        % dependence is set by default grid and there is no time or
        % position dependce.
        dT = obj.per_scat_kern; 
    end
    
    
    function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        error('Function not yet implemented!')
    end

    
    function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        error('Function not yet implemented!')
    end
    
    
    function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        error('Function not yet implemented!')
    end
    
    
    function h_i = getOneParticleEV(obj, charIdx, t, x, rapid)
        % =================================================================
        % Purpose : Returns one-particle eigenvalues of the charIdx'th charge. 
        % Input :   charIdx -- array of charges to calculate
        %           t     -- time (scalar)
        %           x     -- x-coordinate (can be scalar or vector)
        %           rapid -- rapid-coordinate (can be scalar or vector)
        % Output:   h_i   -- i'th order eigenvalue
        % =================================================================
        switch charIdx
        case 0 % eigenvalue of number operator
            h_i = obj.getBareTopologicalCharge();
            h_i = repmat(h_i, [length(rapid), 1]);

        case 1 % eigenvalue of momentum operator
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);

        case 2 % eigenvalue of Hamiltonian
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid);

        otherwise 
            % Higher order charges must be implemented in specific model
            error(['Eigenvalue ' num2str(charIdx) ' not implmented!'])
        end
    end


    function e_eff = calcEffectiveEnergy(obj, w, t, x)
        % Overloaded from superclass (iFluidCore)
        % w is source term

        % For solving TBA equations, periodic boundaries for magnons are
        % necessary
        kernel      = 1/(2*pi)*obj.getPeriodicScatteringRapidDeriv(0,0,0,0,0,0); % arguments dont matter
        eta         = obj.getTypeSigns();
        

        e_eff       = fluidcell(w);
        e_eff_old   = fluidcell(w);
        error_G     = 1;
        count       = 0;
        
        % Solve TBA eq. for epsilon(k) by iteration:
        % Using epsilonk_old, update epsilon_k until convergence is reached 
        % (difference between epsilonk_old and epsilon_k becomes less than tol)
    
        switch obj.TBA_solver

        case 'Picard'
            while any(error_G > obj.tolerance) && count < obj.maxcount % might change this
                
                % calculate epsilon(k) from integral equation using epsilonk_old
                % i.e. update epsilon^[n] via epsilon^[n-1]      
                e_eff   = w + kernel*(obj.rapid_w.*eta.*obj.getFreeEnergy(e_eff)); % free energy has negative sign

                error_G = squeeze(double( sum(abs(e_eff - e_eff_old), [1, 2]) ))/obj.N/obj.M;
                e_eff_old = e_eff;

                count   = count+1;
            end

        case 'Newton'
            while any(error_G > obj.tolerance) && count < obj.maxcount % might change this
                
                I       = fluidcell.eye(obj.N, obj.Ntypes);
                G       = e_eff - ( w + kernel*(obj.rapid_w.*eta.*obj.getFreeEnergy(e_eff)) ); % free energy has negative sign
                dGde    = I - kernel*(obj.rapid_w.*eta.*(I./(1+exp(e_eff)) ));
                e_eff   = e_eff - (dGde+eps)\G;               

                error_G = squeeze(double( sum(abs(G), [1, 2]) ))/obj.N/obj.M;

                count   = count+1;
            end           
        otherwise
            error('The specified TBA solver has not been implemented.')
        end

        fprintf('TBA equation converged with error %d in %d iterations.\n', max(error_G), count);
    end
    
       
    function [f1, f2] = calcFreeEnergyDensity(obj, e_eff, T, rho, rhoS)
        % returns f/T, i.e. free energy in units of temperature
        eta     = obj.getTypeSigns();
        masses  = obj.getMasses(obj.type_grid);
        ebare   = masses.*cosh(obj.rapid_grid);
        mu      = obj.couplings{1,1}(0,obj.x_grid);
        rhoH    = rhoS - rho;
        q       = obj.getOneParticleEV(0, 0, obj.x_grid, obj.rapid_grid);

        f1      = 1/(2*pi)*sum(obj.rapid_w.*sum(eta.*ebare.*obj.getFreeEnergy(e_eff), 3), 1);
        f2      = 1/T*sum(sum(obj.rapid_w .* (rho.*ebare - T*rho.*log(1 + rhoH./(rho+eps)) - T*rhoH.*log(1 + rho./(rhoH+eps)) - rho.*mu.*q) , 1), 3);

        f1      = squeeze(double(f1));
        f2      = squeeze(double(f2));
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
        % overload from superclass (iFluidCore) for numerically more stable
        % version
        
        theta = exp(-e_eff)./( exp(-e_eff) + 1);
    end
    
    
    function Q_dr = applyDressing(obj, Q, theta, t)
        % =================================================================
        % Purpose : Dresses quantity Q by solving system of linear eqs.
        % Input :   Q     -- Quantity to be dressed
        %           theta -- filling function (fluidtensor)
        %           t     -- time (scalar)
        % Output:   Q_dr  -- Dressed quantity (fluidtensor)
        % =================================================================
        if ~isa(Q, 'fluidcell')
            Q = fluidcell(Q);
        end
        
        if size(Q,1) == 1
            X = repmat(double(Q), obj.N, 1, obj.Ntypes);
            Q = fluidcell(X);
        end
        
        % Calculate dressing operator
        eta     = obj.getTypeSigns();
        T       = 1/(2*pi)*obj.getScatteringRapidDeriv(0,0,0,0,0,0); % arguments dont matter
        I       = fluidcell.eye(obj.N, obj.Ntypes);
        U       = I./eta - T.*transpose(obj.rapid_w.*theta);
        
        % We now have the equation Q = U*Q_dr. Therefore we solve for Q_dr
        % using the '\' operation.
        Q_dr     = U\Q;

        % For kernels, the solution to the dressing equation must be transposed 
        if size(Q_dr, 1) == size(Q_dr, 4)
            Q_dr = Q_dr.t();
        end
    end
    
      
end % end public methods

    
end % end classdef