classdef SineGordonModel_reflectionless < iFluidCore
    % iFluid implmentation of TBA functions for sine-Gordon model
    %
    % ### Couplings are {mu}
    %
properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies    = 'fermion'; 
    
    N_bre           = [];  % number of breather types
    M_sol           = 1;   % mass of soliton 
    xi              = []; % renormalized coupling

    scat_kern       = []; % scattering kernel

end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = SineGordonModel_reflectionless(x_grid, rapid_grid, rapid_w, N_bre, couplings, Options)   
        % Set default parameters
        if nargin < 5
            Options = struct;
        end
        
        % The coupling is specified by the number of breather types N_bre
        Ntypes  = N_bre + 2; % plus soliton and anti-soliton 

        % Call superclass constructor
        obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
        
        obj.N_bre = N_bre;
        obj.xi    = 1/( N_bre + 1 );

        obj.scat_kern = fluidcell( obj.getTotalScatteringKernel() );
    end
    
    
    %% Implementation of abstract equations
    
    function M = getMasses(obj, type)
        % get mass of breather/soliton with the given type index
        % type = 1:N_bre --> breathers
        % type = N_bre + 1 --> soliton
        % type = N_bre + 2 -- > anti-soliton
  
        M = 2*obj.M_sol*sin(type*pi*obj.xi/2); % breather masses
        
        M(type == obj.N_bre + 1) = obj.M_sol; % soliton mass
        M(type == obj.N_bre + 2) = obj.M_sol; % anti-soliton mass
    end


    function xi = getXi(obj)
        xi = obj.xi;
    end


    function eta = getTypeSigns(obj)
        % Returns parity of each quasiparticle species
        eta                 = zeros(1,1,obj.Ntypes);
        eta(1:obj.N_bre)    = 1; % breathers
        eta(obj.N_bre + 1)  = 1; % soliton
        eta(obj.N_bre + 2)  = 1; % anti-soliton

        eta                 = fluidcell(eta);
    end


    function qbare = getBareTopologicalCharge(obj)
        qbare               = zeros(1,1,obj.Ntypes); % breathers have no charge
        qbare(obj.N_bre+1)  = 1; % soliton
        qbare(obj.N_bre+2) = -1; % anti-soliton
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
        % i is breather index

        rapid   = obj.rapid_grid - obj.rapid_grid';
        
        Phi     = -4*cos(i*pi*obj.xi/2)*cosh(rapid)./(cos(i*pi*obj.xi) + cosh(2*rapid));

        for k = 1:i-1
            Phi = Phi - 2*cos((i-2*k)*pi*obj.xi/2)./(cosh(rapid) - sin((i-2*k)*pi*obj.xi/2));
        end
    end


    function Phi = getScatteringKernel_BB(obj, i, j)
        % Returns Breather-Breather scattering kernel as N x N matrix
        % i and j are breather index

        rapid = obj.rapid_grid - obj.rapid_grid';

        Phi = 4*sin((i+j)*pi*obj.xi/2)*cosh(rapid)./(cos((i+j)*pi*obj.xi) - cosh(2*rapid) + eps) + ...
              4*sin((i-j)*pi*obj.xi/2)*cosh(rapid)./(cos((i-j)*pi*obj.xi) - cosh(2*rapid) + eps);

        for k = 1:j-1
            Phi = Phi - 2*sin((j-i-2*k)*pi*obj.xi/2)./(cos((j-i-2*k)*pi*obj.xi/2) - cosh(rapid) + eps) - ...
                        2*sin((j+i-2*k)*pi*obj.xi/2)./(cos((j+i-2*k)*pi*obj.xi/2) + cosh(rapid) + eps);
        end
    end


    function Phi_tot = getTotalScatteringKernel(obj)

        Phi_tot = zeros(obj.N, obj.N, obj.Ntypes, obj.Ntypes);

        % calculate soliton-soliton scattering (assume for now that this is
        % given by the self-scattering \Phi_0)
        Phi     = obj.getScatteringKernel_SS();
        Phi_tot(:, :, obj.N_bre+1, obj.N_bre+1) = Phi; % SS
        Phi_tot(:, :, obj.N_bre+1, obj.N_bre+2) = Phi; % SA
        Phi_tot(:, :, obj.N_bre+2, obj.N_bre+1) = Phi; % AS
        Phi_tot(:, :, obj.N_bre+2, obj.N_bre+2) = Phi; % AA

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

            Phi_tot(:, :, obj.N_bre+2, i) = Phi;
            Phi_tot(:, :, i, obj.N_bre+2) = Phi;
        end

        Phi_tot = fluidcell(permute(Phi_tot, [1 5 3 2 4]));

    end


    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % For now, scattering kernel is pre-calculated
        dT = obj.scat_kern; 
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
        kernel      = 1/(2*pi)*obj.getScatteringRapidDeriv(0,0,0,0,0,0); % arguments dont matter
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