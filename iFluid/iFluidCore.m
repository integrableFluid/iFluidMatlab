classdef iFluidCore < handle
    % Superclass for solving GHD dynamical equations.
    % Extends handle class --> properties can be changed with setter
    % methods while retaining GHD object.
    % This class contains abstract methods, which are model specific,
    % whereby they must be implemented in the corresponding subclass.
    % 
    % NOTE: This version employs an auto-derivative feature, whereby only
    % the model parameters (energy, momentum, scattering) along with the
    % couplings must be specified - their derivatives are automatically
    % taken. The model parameters, couplings, and associated derivatives
    % are all stores as anonymous functions in cell-arrays.
    % 
    % In case the auto-derivative fails (feature is not very robust), one
    % must supply coupling-derivatives along with the couplings in a cell
    % array, while the derivatives of model parameters must be implmented 
    % as overloaded class functions.
    %
    % Several options exists for solvers and auto-derivatives. The default
    % options can be found in the properties of this class. To change
    % options one must supply the constructor with a struct containing the
    % options-parameters one wishes to change.
    %
    % This class and any subclass of it must follow the following
    % convention for indeces:
    %   Index 1: Main rapidity index, size N.
    %   Index 2: Secondary rapid index, used for kernels.
    %   Index 3: Main type index, size Ntypes.
    %   Index 4: Secondary type index, used for kernels.
    %   Index 5: Spatial index, size M.
    %
    
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
    
    
    % Cell array with anonymous functions
    couplings       = []; % @(t,x) coupling #i
    
    % Optional parameters (default values specified here). 
    tolerance       = 1e-6;     % Tolerance for TBA solution
    maxcount        = 100;      % Max interations for TBA solution
    homoEvol        = false;

end % end protected properties


properties (Abstract, Access = protected)
    % These properties must be specified in the specific implementation

    % Species of quasiparticle
    quasiSpecies
    
end % end abstract properties

methods (Abstract, Access = public)
    % These methods must be implemented in the extended class
    ebare   = getBareEnergy(obj, t, x, rapid, type)
    pbare   = getBareMomentum(obj, t, x, rapid, type)
    de      = calcEnergyRapidDeriv(obj, t, x, rapid, type)
    dp      = calcMomentumRapidDeriv(obj, t, x, rapid, type)
    dT      = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
    de      = calcEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    dp      = calcMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    dT      = calcScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
    
    
end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options)        
        % Check format of input
        if ~isvector(x_grid) || ~isvector(rapid_grid)
            error('Input has wrong format!')
        end
        if ~iscell( couplings )
            error('couplings must be cell array of anonymous functions!')
        end
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = Ntypes;
        
        % Reshape grids to right format
        obj.x_grid      = reshape(x_grid, 1, 1, 1, 1, obj.M); % 5th index is space
        obj.rapid_grid  = reshape(rapid_grid, obj.N, 1); % 1st index is rapidity
        obj.type_grid   = reshape( 1:Ntypes, 1, 1, Ntypes ); % Types are 3rd index
        obj.rapid_w     = reshape( rapid_w , length(rapid_w), 1); % 1st index is rapidity
        
        % Copy fields of Options struct into class properties
        if nargin > 4
            fn = fieldnames(Options);
            for i = 1:length(fn)                      
                if isprop(obj, fn{i}) % only copy field if defined among properties
                    eval(['obj.',fn{i},' = Options.',fn{i},';']);
                end
            end
        end
        
        % Fill out cell arrays with anonymous functions
        obj.setCouplings(couplings);
        
    end
    
    
    function setCouplings(obj, couplings)
        if ~iscell( couplings )
            error('couplings must be cell array of anonymous functions!')
        end
        if size(couplings, 1) ~= 3
            obj.couplings = couplings(1,:);
        else  
            obj.couplings = couplings;
        end
    end
    
    
    function couplings = getCouplings(obj)
        couplings = obj.couplings;
    end


    function [x_grid, rapid_grid, type_grid, rapid_w] = getGrids(obj)
        x_grid      = obj.x_grid;
        rapid_grid  = obj.rapid_grid;
        type_grid   = obj.type_grid;
        rapid_w     = obj.rapid_w;
    end
    
    
    function f = getStatFactor(obj, theta)
        switch obj.quasiSpecies
        case 'fermion'
            f = 1 - theta;
        case 'boson'
            f = 1 + theta;
        case 'classical'
            f = 1;
        case 'radiative'
            f = theta;
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end
    end

    
    function F = getFreeEnergy(obj, e_eff)
        switch obj.quasiSpecies
        case 'fermion'
            F = -log( 1 + exp(-e_eff));
        case 'boson'
            F = log( 1 - exp(-e_eff));
        case 'classical'
            F = -exp(-e_eff);
        case 'radiative'
            F = log( e_eff );
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end 
    end
    
    
    function h_i = getOneParticleEV(obj, charIdx, t, x, rapid)
        switch charIdx
        case 0 % eigenvalue of number operator
            h_i = repmat(obj.type_grid, length(rapid), 1);
        case 1 % eigenvalue of momentum operator
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of Hamiltonian
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid);
        otherwise 
            error(['Eigenvalue ' num2str(charIdx) ' not implmented!'])
        end
    end
    
       
    function [rho, rhoS] = transform2rho(obj, theta, t_array)
        % Transforms from occupation number basis (theta) to
        % occupied density of states basis (rho).
        % NOTE: if couplings are time-dependent, t = 0 should be used for
        % units to stay fixed.
        
        Nsteps  = length(theta); % number of time steps
        rho     = cell(1,Nsteps);
        rhoS    = cell(1,Nsteps);
        for n = 1:Nsteps % Transform for each time step
            
            if Nsteps == 1
                theta_n = theta;
            else
                theta_n = theta{n};
            end
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            dp_dr   = obj.applyDressing(dp, theta_n, t);
            
            rhoS{n} = 1/(2*pi) * dp_dr;
            rho{n}  = theta_n.*rhoS{n}; 
        end
        
        if Nsteps == 1
            rho     = rho{1};
            rhoS    = rhoS{1};
        end
    end
    
    
    function [theta, rhoS] = transform2theta(obj, rho, t_array)
        % Transforms from occupied density of states basis (rho) to
        % occupation number basis (theta).
        % NOTE: if couplings are time-dependent, t = 0 should be used for
        % units to stay fixed.
        
        Nsteps  = length(rho); % number of time steps
        theta   = cell(1,Nsteps);
        rhoS    = cell(1,Nsteps);
        for n = 1:Nsteps % Transform for each time step
            if Nsteps == 1
                rho_n = rho;
            else
                rho_n = rho{n};
            end
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            kernel  = 1/(2*pi) * obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid); 
            
            rhoS{n} = dp/(2*pi) - kernel*(obj.rapid_w.*rho_n);
            theta{n}= rho_n./rhoS{n};
        end
        
        if Nsteps == 1
            theta   = theta{1};
            rhoS    = rhoS{1};
        end
    end
    

    function [q, j] = calcCharges(obj, theta, c_idx, t_array)
        % Calculate charges, q_i, and associated currents, j_i, where i is
        % entry in vector c_idx.
        
        if iscell(theta)
            Nsteps = length(theta); % number of time steps
        else
            Nsteps = 1;
        end
        
        Ncharg  = length(c_idx); 
        q       = zeros(obj.M, Nsteps, Ncharg);
        j       = zeros(obj.M, Nsteps, Ncharg);
    
        for i = 1:Nsteps
            if Nsteps == 1
                theta_i = theta;
            else
                theta_i = theta{i};
            end
            
            if nargin < 4
                t = 0;
            else
                t = t_array(i);
            end
            
            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            dE      = obj.calcEnergyRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);

            for n = 1:Ncharg
                hn          = obj.getOneParticleEV( c_idx(n), t, obj.x_grid, obj.rapid_grid);               
                hn_dr       = obj.applyDressing(hn, theta_i, t);
                
                q(:,i,n)    = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*theta_i.*hn_dr) , 3) , 1));
                j(:,i,n)    = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dE.*theta_i.*hn_dr) , 3) , 1)); 
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % Calculate initial filling function, given by the equilibrium TBA
        % for some temperature, T, and couplings.
        % If no couplings are passed to method, use already set couplings.
        
        if nargin == 3 % use input TBA_couplings
            % Save old couplings 
            couplings_old   = obj.getCouplings();
            
            % Set new couplings
            obj.setCouplings(TBA_couplings);
        end
            
        e_eff = obj.calcEffectiveEnergy(T, 0, obj.x_grid, obj.rapid_grid);
        theta = obj.calcFillingFraction(e_eff);
        
        if nargin == 3
            % Return to old couplings
            obj.setCouplings(couplings_old);
        end
    end


    function [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid, type)        
        % Calculates velocities
        de_dr   = obj.applyDressing(obj.calcEnergyRapidDeriv(t, x, rapid, type), theta, t);
        dp_dr   = obj.applyDressing(obj.calcMomentumRapidDeriv(t, x, rapid, type), theta, t);
        
        v_eff   = de_dr./dp_dr;
         
        % Calculate acceleration from inhomogenous couplings
        a_eff   = 0;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = iFluidTensor( zeros(size( v_eff )) );
            return
        end
        
        % Calculate contribution for each coupling
        for coupIdx = 1:size(obj.couplings,2)            
            dT      = obj.calcScatteringCouplingDeriv(coupIdx, t, x, rapid, obj.rapid_grid, type, obj.type_grid);
            accKern = 1/(2*pi) * dT.*transpose(obj.rapid_w .* theta);
            
            % if time deriv of coupling exist compute f
            if ~isempty(obj.couplings{2,coupIdx}) 
                f       = -obj.calcMomentumCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*dp_dr;
                f_dr    = obj.applyDressing(f, theta, t);
                a_eff   = a_eff + obj.couplings{2,coupIdx}(t,x).*f_dr;
            end
            
            % if spacial deriv of coupling exist compute Lambda
            if ~isempty(obj.couplings{3,coupIdx}) 
                L       = -obj.calcEnergyCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*de_dr;
                L_dr    = obj.applyDressing(L, theta, t);
                a_eff   = a_eff + obj.couplings{3,coupIdx}(t,x).*L_dr;
            end
        end
        
     
        a_eff   = a_eff./dp_dr;
    end
    
    function Q_dr = applyDressing(obj, Q, theta, t)
        % Applies dressing to quantity Q. Q must have first index of
        % dimension N!
        
        if ~isa(Q, 'iFluidTensor')
            Q = iFluidTensor(Q);
        end
        
        if size(Q,1) == 1
            X = repmat(double(Q), obj.N, 1, obj.Ntypes);
            Q = iFluidTensor(X);
        end
        
        % Calculate dressing operator
        kernel  = 1/(2*pi)*obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid);
        
        I_rapid = eye(obj.N);
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        identity= I_rapid.*I_type;

        U       = identity + kernel.*transpose(obj.rapid_w.*theta); 
        
        % We now have the equation Q = U*Q_dr. Therefore we solve for Q_dr
        % using the '\' operation.
        Q_dr     = U\Q;
           
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
        switch obj.quasiSpecies
        case 'fermion'
            theta = 1./( exp(e_eff) + 1);
        case 'boson'
            theta = 1./( exp(e_eff) - 1);
        case 'classical'
            theta = exp(-e_eff);
        case 'radiative'
            theta = 1./exp(e_eff);
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end
    end
    
    
    function e_eff = calcEffectiveEnergy(obj, T, t, x, rapid)
        % Calculates the dressed energy per particle epsilon(k), from which the
        % preasure and later the filling factor theta can be derived.
        % This is achieved by iteratively solving the TBA equation.
        ebare       = obj.getBareEnergy(t, x, obj.rapid_grid, obj.type_grid); 
        kernel      = 1/(2*pi)*obj.calcScatteringRapidDeriv(t, x, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid );
        
        if isa(T, 'function_handle')
            T = T(x);
        end
        
        e_eff       = iFluidTensor(obj.N, 1, obj.Ntypes, 1, obj.M);
        e_eff_old   = iFluidTensor(obj.N, 1, obj.Ntypes, 1, obj.M);
        error_rel   = 1;
        count       = 0;
        
        % Solve TBA eq. for epsilon(k) by iteration:
        % Using epsilonk_old, update epsilon_k until convergence is reached 
        % (difference between epsilonk_old and epsilon_k becomes less than tol)
        while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
            
            % calculate epsilon(k) from integral equation using epsilonk_old
            % i.e. update epsilon^[n] via epsilon^[n-1]            
            e_eff       = ebare./T - kernel*(obj.rapid_w .* obj.getFreeEnergy(e_eff_old));
            
            % calculate error
            v1          = flatten(e_eff);
            v2          = flatten(e_eff_old);
            
            sumeff      = sum( v1.^2 ,1);            
            error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
            e_eff_old   = e_eff;

            count       = count+1;
        end
    end
    
    
    
end % end public methods


end % end classdef