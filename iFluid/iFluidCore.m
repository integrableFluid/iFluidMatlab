classdef iFluidCore < handle
    % Superclass for solving general equations of the thermodynamic Bethe
    % ansatz.
    % This class contains abstract methods, which are model specific,
    % whereby they must be implemented in subclasses.
    % 
    % =========== How to extend this class ========
    %
    % ** Step 1 ** Create new model 
    %   
    %   myModel < iFluidCore
    %
    %
    % ** Step 2 ** Implement the abstract properties and methods 
    %   
    %   quasiSpecies    % can be either 'fermion', 'boson', 'classical'
    %
    %   getBareEnergy(obj, t, x, rapid, type)
    %   getBareMomentum(obj, t, x, rapid, type)
    %   getEnergyRapidDeriv(obj, t, x, rapid, type)
    %   getMomentumRapidDeriv(obj, t, x, rapid, type)
    %   getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
    %   getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    %   getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    %   getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
    %
    %
    % ** Step 3 ** Write constructor calling the super-constructor
    %
    %   function obj = myModel(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options)
    %       obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
    %   end
    %
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
    
    rapid_aux       = []; % auxillary grids
    type_aux        = [];

    
    % Cell array with anonymous functions
    couplings       = []; % @(t,x) coupling #i
    
    % Optional parameters (default values specified here). 
    tolerance       = 1e-10;     % Tolerance for TBA solution
    maxcount        = 1000;     % Max interations for TBA solution
    TBA_solver      = 'Newton'; % Scheme for solving implicit TBA eq
    homoEvol        = false;    % Homogeneous evolution. if true, dont calculate effective acceleration

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
    de      = getEnergyRapidDeriv(obj, t, x, rapid, type)
    dp      = getMomentumRapidDeriv(obj, t, x, rapid, type)
    dT      = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
    de      = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    dp      = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
    dT      = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
    
    
end % end protected abstract methods


methods (Access = public)
    
    function obj = iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options)        
        % =================================================================
        % Purpose : Superclass constructor. 
        % Input :   x_grid     -- Vector containing position gridpoints
        %           rapid_grid -- Vector containing rapidity gridpoints
        %           rapid_w    -- Rapidity quadrature weights for integrals
        %           couplings  -- Cell array of anonymous functions for the
        %                          couplings of the specific function. 
        %                          May also include spatial and temporal
        %                          derivatives of the couplings.
        %           Ntypes     -- Number of quasiparticle species/types
        %           Options    -- Struct containing values for optional 
        %                          parameters (see class Properties) 
        % Output:   obj        -- iFluidCore object
        % ================================================================= 
        
        % Check format of input
        if ~isvector(x_grid) || ~isvector(rapid_grid)
            error('Input has wrong format!')
        end
        if ~iscell( couplings )
            error('couplings must be cell array of anonymous functions!')
        end
        
        % Reshape and store grids
        obj.setGrids(x_grid, rapid_grid, rapid_w, Ntypes);

        % Fill out cell arrays with anonymous functions
        obj.setCouplings(couplings);
        
        % Copy fields of Options struct into class properties
        if nargin > 4
            fn = fieldnames(Options);
            for i = 1:length(fn)                      
                if isprop(obj, fn{i}) % only copy field if defined among properties
                    eval(['obj.',fn{i},' = Options.',fn{i},';']);
                end
            end
        end
        
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
    
    
    function setGrids(obj, x_grid, rapid_grid, rapid_w, Ntypes)
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = Ntypes;
        
        
        % Reshape grids to right format
        obj.rapid_grid  = reshape(rapid_grid, obj.N, 1); % 1st index is rapidity
        obj.x_grid      = reshape(x_grid, 1, obj.M); % 2nd index is space
        obj.type_grid   = reshape( 1:Ntypes, 1, 1, Ntypes ); % Types are 3rd index
        
        obj.rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        obj.type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        
        obj.rapid_w     = reshape( rapid_w , length(rapid_w), 1); % 1st index is rapidity
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

    
    function [F, dF_de] = getFreeEnergy(obj, e_eff)
        switch obj.quasiSpecies
        case 'fermion'
            F       = -log( 1 + exp(-e_eff));
            dF_de   = (exp(e_eff) + 1).^(-1);
        case 'boson'
            F       = log( 1 - exp(-e_eff));
            dF_de   = (exp(e_eff) - 1).^(-1);
        case 'classical'
            F       = -exp(-e_eff);
            dF_de   = exp(-e_eff);
        case 'radiative'
            F       = log( e_eff );
            dF_de   = e_eff.^(-1);
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end 
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
            h_i = ones(length(rapid), 1, obj.Ntypes);
        case 1 % eigenvalue of momentum operator
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of Hamiltonian
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid);
        otherwise 
            % Higher order charges must be implemented in specific model
            error(['Eigenvalue ' num2str(charIdx) ' not implmented!'])
        end
    end
    
       
    function [rho, rhoS] = transform2rho(obj, theta, t_array)
        % =================================================================
        % Purpose : Transforms filling function into root density.
        % Input :   theta   -- filling function (cell array of ifluidcell)
        %           t_array -- array of times corresponding to theta
        % Output:   rho     -- root density (cell array of ifluidcell)
        %           rhoS    -- density of states (cell array of ifluidcell)
        % =================================================================   
        
        if iscell(theta)
            Nsteps = length(theta); % number of time steps
        else
            Nsteps = 1;
        end
        
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

            dp      = obj.getMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
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
        % =================================================================
        % Purpose : Transforms root density into filling function.
        % Input :   rho     -- root density (cell array of ifluidcell)
        %           t_array -- array of times corresponding to theta
        % Output:   theta   -- filling function (cell array of ifluidcell)
        %           rhoS    -- density of states (cell array of ifluidcell)
        % =================================================================
        
        if iscell(rho)
            Nsteps = length(rho); % number of time steps
        else
            Nsteps = 1;
        end
        
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

            dp      = obj.getMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            kernel  = 1/(2*pi) * obj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux); 
            
            rhoS{n} = dp/(2*pi) - kernel*(obj.rapid_w.*rho_n);
            theta{n}= rho_n./rhoS{n};
        end
        
        if Nsteps == 1
            theta   = theta{1};
            rhoS    = rhoS{1};
        end
    end
    

    function [q, j, Vq, Vj] = calcCharges(obj, c_idx, theta, t_array, varargin)
        % =================================================================
        % Purpose : Calculates expectation values of charge densities and
        %           associated currents.
        % Input :   c_idx   -- charge indices
        %           theta   -- filling function (cell array of ifluidcell)
        %           t_array -- array of times corresponding to theta
        %           varargin-- optional arguments:
        %                       calc_formfactors (default = false)
        %                       sum_types (default = true)
        % Output:   q       -- charge exp. vals for each t in t_array
        %           j       -- current exp. vals for each t in t_array
        %           Vq      -- one-particle charge form factors
        %           Vq      -- one-particle current form factors
        % ================================================================= 
        
        parser = inputParser;

        addParameter(parser, 'calc_formfactors', false, @(x) islogical(x));
        addParameter(parser, 'sum_types', true, @(x) islogical(x));

        parse(parser, varargin{:});
        
        
        Ncharg  = length(c_idx); 
        Nsteps  = length(t_array);
        
        if parser.Results.sum_types
            q       = zeros(obj.M, Nsteps, Ncharg);
            j       = zeros(obj.M, Nsteps, Ncharg);
            Vq      = cell(Nsteps, Ncharg);
            Vj      = cell(Nsteps, Ncharg);
        else
            q       = zeros(obj.M, Nsteps, Ncharg, obj.Ntypes);
            j       = zeros(obj.M, Nsteps, Ncharg, obj.Ntypes);
            Vq      = cell(Nsteps, Ncharg);
            Vj      = cell(Nsteps, Ncharg);
        end
    
        for i = 1:Nsteps
            if Nsteps == 1
                theta_i = theta;
            else
                theta_i = theta{i};
            end
            
            t       = t_array(i);
            
            dp      = obj.getMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            dE      = obj.getEnergyRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            
            if parser.Results.calc_formfactors
                v_eff = obj.calcEffectiveVelocities(theta_i, t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            end

            for n = 1:Ncharg
                hn          = obj.getOneParticleEV( c_idx(n), t, obj.x_grid, obj.rapid_grid);               
                hn_dr       = obj.applyDressing(hn, theta_i, t);
                
                if parser.Results.sum_types
                q(:,i,n)    = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*theta_i.*hn_dr) , 3) , 1));
                j(:,i,n)    = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dE.*theta_i.*hn_dr) , 3) , 1));
                else
                q(:,i,n,:)  = 1/(2*pi) * squeeze(sum( obj.rapid_w .* double(dp.*theta_i.*hn_dr), 1));
                j(:,i,n,:)  = 1/(2*pi) * squeeze(sum( obj.rapid_w .* double(dE.*theta_i.*hn_dr), 1));
                end
                if parser.Results.calc_formfactors
                    Vq{i,n} = hn_dr;
                    Vj{i,n} = v_eff.*hn_dr;
                end
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % =================================================================
        % Purpose : Calculates thermal state for the given couplings and
        %           temperature.
        % Input :   T     -- Temperature (scalar or @(x))
        %           TBA_couplings -- (optional) cell array of couplings. If
        %                            none are passed, use already set
        %                            couplings of model.
        % Output:   theta -- Filling function of thermal state 
        %           e_eff -- Pesudo-energy of thermal state 
        % ================================================================= 
        if nargin == 3 % use input TBA_couplings
            % Save old couplings 
            couplings_old   = obj.getCouplings();
            
            % Set new couplings
            obj.setCouplings(TBA_couplings);
        end
            
        ebare   = obj.getBareEnergy(0, obj.x_grid, obj.rapid_grid, obj.type_grid);
        if isa(T, 'function_handle')
            T       = T(obj.x_grid);
        end
        
        w       = ebare./T;
        
        e_eff   = obj.calcEffectiveEnergy(w, 0, obj.x_grid);
        theta   = obj.calcFillingFraction(e_eff);
        
        if nargin == 3
            % Return to old couplings
            obj.setCouplings(couplings_old);
        end
    end

    
    function [v_eff, a_eff, de_dr, dp_dr] = calcEffectiveVelocities(obj, theta, t, varargin)        
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration of
        %           quasiparticles.
        % Input :   theta -- filling function (ifluidcell)
        %           t     -- time (scalar)
        % Output:   v_eff -- Effective velocity (ifluidcell)
        %           a_eff -- Effective acceleration (ifluidcell)
        % =================================================================
        if length(varargin) == 0
            % No extra arguments given, use default grids
            x       = obj.x_grid;
            rapid   = obj.rapid_grid;
            type    = obj.type_grid;
            
            [v_eff, a_eff, de_dr, dp_dr] = obj.calcVelocitiesNormal(theta, t, x, rapid, type);  
            
        elseif length(varargin) == 3
            % Grids given as extra arguments
            x       = varargin{1};
            rapid   = varargin{2};
            type    = varargin{3};
            
            [v_eff, a_eff, de_dr, dp_dr] = obj.calcVelocitiesNormal(theta, t, x, rapid, type);  
            
        elseif length(varargin) == 1
            % Dressing kernel given as extra argument. Calculate dressing
            % using this kernel instead of '\' operation.
            D       = varargin{1};
            
            [v_eff, a_eff, de_dr, dp_dr] = obj.calcVelocitiesFast(theta, t, D);  
        end
    end
    

    function [v_eff, a_eff, de_dr, dp_dr] = calcVelocitiesNormal(obj, theta, t, x, rapid, type)        
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration of
        %           quasiparticles.
        % Input :   theta -- filling function (fluidcell)
        %           t     -- time (scalar)
        %           x     -- x-coordinate (can be scalar or vector)
        %           rapid -- rapid-coordinate (can be scalar or vector)
        %           type  -- quasiparticle type (can be scalar or vector)
        % Output:   v_eff -- Effective velocity (fluidcell)
        %           a_eff -- Effective acceleration (fluidcell)
        % =================================================================
        
        de_dr   = obj.applyDressing(obj.getEnergyRapidDeriv(t, x, rapid, type), theta, t);
        dp_dr   = obj.applyDressing(obj.getMomentumRapidDeriv(t, x, rapid, type), theta, t);
        
        v_eff   = de_dr./dp_dr;
         
        % Calculate acceleration from inhomogenous couplings
        a_eff   = 0;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = fluidcell.zeros( size(v_eff) );
            return
        end
        
        % Calculate contribution for each coupling
        for coupIdx = 1:size(obj.couplings,2)            
            dT      = obj.getScatteringCouplingDeriv(coupIdx, t, x, rapid, obj.rapid_aux, type, obj.type_aux);
            accKern = 1/(2*pi) * dT.*transpose(obj.rapid_w .* theta);
            
            % if time deriv of coupling exist, compute f
            if ~isempty(obj.couplings{2,coupIdx}) 
                f       = -obj.getMomentumCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*dp_dr;
                f_dr    = obj.applyDressing(f, theta, t);
                a_eff   = a_eff + obj.couplings{2,coupIdx}(t,x).*f_dr;
            end
            
            % if spacial deriv of coupling exist, compute Lambda
            if ~isempty(obj.couplings{3,coupIdx}) 
                L       = -obj.getEnergyCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*de_dr;
                L_dr    = obj.applyDressing(L, theta, t);
                a_eff   = a_eff + obj.couplings{3,coupIdx}(t,x).*L_dr;
            end
        end
        
     
        a_eff   = a_eff./dp_dr;
    end
    
    
    function [v_eff, a_eff, de_dr, dp_dr] = calcVelocitiesFast(obj, theta, t, D)        
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration of
        %           quasiparticles.
        % Input :   theta -- filling function (fluidcell)
        %           t     -- time (scalar)
        %           D     -- dressing operator (fluidcell)
        % Output:   v_eff -- Effective velocity (fluidcell)
        %           a_eff -- Effective acceleration (fluidcell)
        % =================================================================
        x       = obj.x_grid;
        rapid   = obj.rapid_grid;
        type    = obj.type_grid;
        
        de      = obj.getEnergyRapidDeriv(t, x, rapid, type);
        dp      = obj.getMomentumRapidDeriv(t, x, rapid, type);
        
        de_dr   = obj.dress(de, D);
        dp_dr   = obj.dress(dp, D);
        
        v_eff   = de_dr./dp_dr;
         
        % Calculate acceleration from inhomogenous couplings
        a_eff   = 0;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = fluidcell.zeros( size(v_eff) );
            return
        end
        
        % Calculate contribution for each coupling
        for coupIdx = 1:size(obj.couplings,2)            
            dT      = obj.getScatteringCouplingDeriv(coupIdx, t, x, rapid, obj.rapid_aux, type, obj.type_aux);
            accKern = 1/(2*pi) * dT.*transpose(obj.rapid_w .* theta);
            
            % if time deriv of coupling exist, compute f
            if ~isempty(obj.couplings{2,coupIdx}) 
                f       = -obj.getMomentumCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*dp_dr;
                f_dr    = obj.dress(f, D);
                a_eff   = a_eff + obj.couplings{2,coupIdx}(t,x).*f_dr;
            end
            
            % if spacial deriv of coupling exist, compute Lambda
            if ~isempty(obj.couplings{3,coupIdx}) 
                L       = -obj.getEnergyCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*de_dr;
                L_dr    = obj.dress(L, D);
                a_eff   = a_eff + obj.couplings{3,coupIdx}(t,x).*L_dr;
            end
        end
        
     
        a_eff   = a_eff./dp_dr;
    end
    
    
    function Q_dr = applyDressing(obj, Q, theta, t)
        % =================================================================
        % Purpose : Dresses quantity Q by solving system of linear eqs.
        % Input :   Q     -- Quantity to be dressed
        %           theta -- filling function (fluidcell)
        %           t     -- time (scalar)
        % Output:   Q_dr  -- Dressed quantity (fluidcell)
        % =================================================================
        if ~isa(Q, 'fluidcell')
            Q = fluidcell(Q);
        end
        
        if size(Q,1) == 1
            X = repmat(double(Q), obj.N, 1, obj.Ntypes);
            Q = fluidcell(X);
        end
        
        % Calculate dressing operator
        kernel  = 1/(2*pi)*obj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux);
        
        I       = fluidcell.eye(obj.N, obj.Ntypes);

        U       = I + kernel.*transpose(obj.rapid_w.*theta);
        
        % We now have the equation Q = U*Q_dr. Therefore we solve for Q_dr
        % using the '\' operation.
        Q_dr     = U\Q;
    end
    
    
    function Q_dr = dress(obj, Q, Uinv)
        % =================================================================
        % Purpose : Dresses a quantity using an inverted dressing kernel.
        %           This is faster when dressing multiple quantities with
        %           the same kernel
        % Input :   Q     -- bare/undressed quantity
        %           Uinv  -- inverted dressing operator (fluidcell)
        % Output:   Q_dr  -- dresse quantity
        % =================================================================

        if ~isa(Q, 'fluidcell')
            Q = fluidcell(Q);
        end
        
        if size(Q,1) == 1
            X = repmat(double(Q), obj.N, 1, obj.Ntypes);
            Q = fluidcell(X);
        end
        
        Q_dr = Uinv*Q;
        
    end
    
    
    function Uinv = calcDressingOperator(obj, theta, t)
        % =================================================================
        % Purpose : Returns the dressing operator used to dress functions.
        % Input :   theta -- filling function (fluidcell)
        %           t     -- time (scalar)
        % Output:   Uinv  -- Dressing operator (fluidcell)
        % =================================================================
        
        kernel  = 1/(2*pi)*obj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux);        
        I       = fluidcell.eye(obj.N, obj.Ntypes);
        U       = I + kernel.*transpose(obj.rapid_w.*theta);
        Uinv    = inv(U);
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
        % =================================================================
        % Purpose : Calculate filling function for given pseudo-energy.
        % Input :   e_eff -- pseudo-energy solution to Yang-Yang equation
        % Output:   theta -- filling function
        % =================================================================

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
    

    function e_eff = calcEffectiveEnergy(obj, w, t, x)
        % =================================================================
        % Purpose : Solves Yang-Yang equation to obtain pseudo-energy of 
        %           local equilibrium state.
        % Input :   w     -- source term corresponding to GGE or thermal 
        %           t     -- time (scalar)
        %           x     -- x-coordinate (can be scalar or vector)
        % Output:   e_eff -- Pesudo-energy corresponding to source term
        % =================================================================

        kernel      = 1/(2*pi)*obj.getScatteringRapidDeriv(t, x, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux );       

        e_eff       = fluidcell(w);
        e_eff_old   = fluidcell(w);
        error_rel   = 1;
        count       = 0;
        
        % Solve Yang-Yang equation by iteration, using either Picard or
        % Newtons method.    
        switch obj.TBA_solver

        case 'Picard'
            while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
                
                % calculate epsilon(k) from integral equation using epsilonk_old
                % i.e. update epsilon^[n] via epsilon^[n-1]      
                e_eff       = w - kernel*(obj.rapid_w .* obj.getFreeEnergy(e_eff_old));
     
                % calculate error
                v1          = flatten(e_eff);
                v2          = flatten(e_eff_old);
    
                sumeff      = sum( v1.^2 ,1);            
                error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
                e_eff_old   = e_eff;
                
                count       = count+1;
            end

        case 'Newton'
            while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
                
                I       = fluidcell.eye(obj.N, obj.Ntypes);
                [F, dF] = obj.getFreeEnergy(e_eff);
    
                G       = e_eff - ( w - kernel*(obj.rapid_w.*F) );
                dGde    = I + kernel*(obj.rapid_w.*(I.*dF) );
    
                e_eff   = e_eff - (dGde+eps)\G;               

     
                % calculate error
                v1          = flatten(e_eff);
                v2          = flatten(e_eff_old);
    
                sumeff      = sum( v1.^2 ,1);            
                error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
                e_eff_old   = e_eff;
                
                count       = count+1;
            end           
        otherwise
            error('The specified TBA solver has not been implemented.')
        end

        fprintf('TBA equation converged with tolerance %d in %d iterations.\n', obj.tolerance, count);
    end
    
    
    function s = calcEntropyDensity(obj, theta, t_array)
       
        [rhoP, rhoS] = obj.transform2rho(theta, t_array);
        
        Nsteps = length(t_array);
        s       = zeros(obj.M, Nsteps);
        
        for i = 1:Nsteps
            s(:,i) = - squeeze(double(sum(sum(obj.rapid_w.*rhoS{i}.*( theta{i}.*log(theta{i}) + (1-theta{i}).*log(1-theta{i}) ) ,1,'omitnan'),3,'omitnan')));
        end

    end
    

    function D = calcDrudeWeight(obj, c_idx, theta, t)
        % =================================================================
        % Purpose : Calculates Drude weight of specified charges.
        % Input :   c_idx   -- charge indices
        %           theta   -- filling function (ifluidcell)
        %           t       -- time (default value 0)
        % Output:   D       -- matrix of Drude weights
        % ================================================================= 
        
        if nargin < 4
            t = 0;
        end

        D      = zeros(length(c_idx), length(c_idx), size(theta, 2));
    
        % calculate common part of Drude weight for all charges
        rho     = obj.transform2rho(theta, t);
        veff    = obj.calcEffectiveVelocities(theta, t);
        f       = obj.getStatFactor(theta);
        D_base  = obj.rapid_w.*rho.*f.*veff.^2;

        % calculate dressed single-particle eigenvalues
        hn_dr   = cell(1, length(c_idx));
        for n = 1:length(c_idx)
            hn      = obj.getOneParticleEV( c_idx(n), t, obj.x_grid, obj.rapid_grid);               
            hn_dr{n}= obj.applyDressing(hn, theta, t);
        end

        % for each combination of charges, calculate Drude weight
        for i = 1:length(c_idx)
        for j = i:length(c_idx)
            D_temp     = sum(sum(D_base.*hn_dr{i}.*hn_dr{j}, 1), 3);
            D(i,j)     = permute(double(D_temp), [1 3 2]);
            D(j,i)     = permute(double(D_temp), [1 3 2]);
        end
        end


    end


    function L = calcOnsagerMatrix(obj, c_idx, theta, t)
        % =================================================================
        % Purpose : Calculates Onsager matrix of specified charges.
        % Input :   c_idx   -- charge indices
        %           theta   -- filling function (ifluidcell)
        %           t       -- time (default value 0)
        % Output:   L       -- Onsager matrix
        % ================================================================= 
        
        if nargin < 4
            t = 0;
        end

        L           = zeros(length(c_idx), length(c_idx), size(theta, 2));
    
        % calculate common part of Drude weight for all charges
        [rho, rhoS] = obj.transform2rho(theta, t);
        veff        = obj.calcEffectiveVelocities(theta, t);
        f           = obj.getStatFactor(theta);
        T           = 1/(2*pi)*obj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux);
        T_dr        = obj.applyDressing(T, theta, t);

        L_base      = 0.5*obj.rapid_w.*permute(obj.rapid_w,[4 2 3 1]).*rho.*f.*rho.t().*f.t().*abs(veff - veff.t()).*T_dr.^2;

        % calculate dressed single-particle eigenvalues
        hn_dr       = cell(1, length(c_idx));
        for n = 1:length(c_idx)
            hn          = obj.getOneParticleEV( c_idx(n), t, obj.x_grid, obj.rapid_grid);               
            hn_dr{n}    = obj.applyDressing(hn, theta, t);
        end

        % for each combination of charges, calculate Drude weight
        for i = 1:length(c_idx)
        for j = i:length(c_idx)
            L_temp      = L_base .* ( hn_dr{i}.t()./rhoS.t() - hn_dr{i}./rhoS ) .* ( hn_dr{j}.t()./rhoS.t() - hn_dr{j}./rhoS );
            L_temp      = sum( L_temp, [1, 3, 4, 5]); % sum over both rapidity and particle indices
            L(i,j)     = permute(double(L_temp), [1 3 2]);
            L(j,i)     = permute(double(L_temp), [1 3 2]);
        end
        end

    end
    
    
end % end public methods


end % end classdef