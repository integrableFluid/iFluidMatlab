classdef LiebLinigerModel < iFluidCore
    % iFluid implmentation of TBA functions for Lieb-Liniger model
    %
    % ### Couplings are { mu, c } ###
    %
    % The units used in this parameterization are as follows:
    % m         = 1/2
    % hbar      = 1
    % g_1d      = 1
    % Lg        = hbar^2/(m*g_1d) = 1 (unit of length)
    % Eg        = 0.5*m*g_1d^2/hbar^2 = 1 (unit of energy)
    % rapidity  = k, whereby p = hbar*k = rapid
    %  
    
properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = LiebLinigerModel(x_grid, rapid_grid, rapid_w, couplings, Options)   
        if nargin < 5
            % If no options, pass empty struct to superclass constructor
            Options = struct;
        end
        
        % Lieb-Liniger model has 1 species of quasi-particles
        Ntypes = 1;
   
        % Call superclass constructor
        obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
    end
    
    
    %% Implementation of abstract equations
    function ebare = getBareEnergy(obj, t, x, rapid, type)
        % First coupling is chemical potential
        ebare = rapid.^2 - obj.couplings{1,1}(t,x);
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        pbare = rapid;
    end
    
    
    function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
        de = 2*rapid;
    end

    
    function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
        dp = repmat(1, length(rapid), 1);
    end
    
    
    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        dT      = -2*obj.couplings{1,2}(t,x)./( (rapid1-rapid2).^2 + obj.couplings{1,2}(t,x).^2 ); 
        
        dT(isnan(dT)) = 0; % removes any NaN
        
        dT      = fluidcell(dT); % Converts to iFluidTensor
    end
    
    
    function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            de = repmat(-1, length(rapid), 1);
        else
            de = 0;
        end
    end

    
    function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
       dp = 0;
    end
    
    
    function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        if coupIdx == 2
            dT = 2*(rapid1-rapid2)./( (rapid1-rapid2).^2 + obj.couplings{1,2}(t,x).^2 );
        else
            dT = 0;
        end
        
        dT(isnan(dT)) = 0; % removes any NaN
        dT = fluidcell(dT); % Converts to iFluidTensor
    end
    
    
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, mu0_guess, setCouplingFlag)
        % =================================================================
        % Purpose : Fits a central chemical potential such that the 
        %           integrated root density of the thermal state produces 
        %           the specified number of atoms.
        % Input :   T       -- Temeprature of system (scalar or @(x))
        %           V_ext   -- External potential as @(x)
        %           Natoms  -- Number of atoms
        %           mu0_guess -- Initial guess for central chemical pot.
        %           setCouplingFlag -- (optional)if true, set coupling to 
        %                               fitted mu
        % Output:   mu0_fit -- Fitted central chemical potential
        % =================================================================
        
        if nargin < 6
            setCouplingFlag = false;
        end
        
        if isempty(V_ext)
            V_ext = obj.couplings{1};
        end
        
        % Fit mu0 to Natoms
        fitfunc     = @(mu0) abs( Natoms - calcNA(obj, mu0, T, V_ext) );
        options     = optimset('Display','iter');
        mu0_fit     = fminsearch(fitfunc, mu0_guess,options);
        
        if setCouplingFlag % adjust couplings to result
            obj.couplings{1,1} = @(t,x) mu0_fit - V_ext(t,x);
        end
        
        function Natoms_fit = calcNA(obj, mu0, T, V_ext)
            % Calculates number of atoms in stationary TBA state given by
            % specied paramters.
            mu_old      = obj.couplings{1,1};
            mu_fit      = @(t,x) mu0 - V_ext(t,x);
            obj.couplings{1,1} = mu_fit;
            
            ebare   = obj.getBareEnergy(0, obj.x_grid, obj.rapid_grid, obj.type_grid);
            if isa(T, 'function_handle')
                T       = T(obj.x_grid);
            end

            w       = ebare./T;
            e_eff   = obj.calcEffectiveEnergy(w, 0, obj.x_grid);
            fill = obj.calcFillingFraction(e_eff);

            obj.couplings{1,1} = mu_old;
            
            dp      = obj.getMomentumRapidDeriv(0, obj.x_grid, obj.rapid_grid, obj.type_grid);

            h0          = ones(obj.N, 1);            
            h0_dr       = obj.applyDressing(h0, fill, 0);
                
            density    = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*fill.*h0_dr) , 3) , 1));
            
            Natoms_fit      = trapz(obj.x_grid, density);
        end % end nested function
    end
    
    
    function [mu0_fit, nu] = fitAtomnumber3D(obj, T, V_ext, Natoms, mu0_guess, N_levels, mu_level, setCouplingFlag)
        % =================================================================
        % Purpose : Estimates the fraction of atoms in transverse excited
        %           states, where each level is treated as a separate
        %           Lieb-Liniger model with chemical potential offset by
        %           the transverse level spacing.
        % Input :   T       -- Temeprature of system (scalar or @(x)).
        %           V_ext   -- External longitudinal potential as @(x).
        %           Natoms  -- Number of atoms.
        %           mu0_guess -- Initial guess for central chemical pot.
        %           N_levels -- Number of transverse levels accounted for.
        %           setCouplingFlag -- (optional)if true, set coupling to 
        %                               fitted mu
        % Output:   mu0_fit -- Fitted central chemical potential.
        %           nu      -- Fitted transverse populations.
        % =================================================================
        
        if nargin < 6
            setCouplingFlag = false;
        end
        
        if isempty(V_ext)
            V_ext = obj.couplings{1};
        end
        
        % Fit mu0 to Natoms
        fitfunc     = @(mu0) abs( Natoms - calcNA(mu0, T, V_ext, N_levels, mu_level) );
        options     = optimset('Display','iter');
        mu0_fit     = fminsearch(fitfunc, mu0_guess,options);
        
        [~, nu]     = calcNA( mu0_fit, T, V_ext, N_levels, mu_level);
        
        if setCouplingFlag % adjust couplings to result
            obj.couplings{1,1} = @(t,x) mu0_fit - V_ext(t,x);
        end
        
        function [Natoms, nu] = calcNA(mu0, T, V_ext, N_levels, mu_level)
            % Calculate total number of atoms
            mu_old  = obj.couplings{1,1};
            NA_i    = zeros(1, N_levels);
            rho_i   = [];
            
            % calculate contribution from each transverse level
            for i = 1:N_levels
                mu_fit      = @(t,x) mu0 - mu_level*(i-1) - V_ext(t,x);
                obj.couplings{1,1} = mu_fit;

                ebare   = obj.getBareEnergy(0, obj.x_grid, obj.rapid_grid, obj.type_grid);
                if isa(T, 'function_handle')
                    T       = T(obj.x_grid);
                end

                w       = ebare./T;
                e_eff   = obj.calcEffectiveEnergy(w, 0, obj.x_grid);
                rho     = calcRho( e_eff );
                density = sum( obj.rapid_w.*rho, 1 );

                NA_i(i) = trapz(obj.x_grid, double(density));
            end

            Natoms  = sum(NA_i); % sum the contributions
            nu      = NA_i/Natoms;
            
            obj.couplings{1,1} = mu_old;
        end % end nested function
        
        function rho = calcRho(e_eff)
            % Calculate the quasiparticle density
            kernel      = obj.getScatteringRapidDeriv(0, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux );
                
            rho         = fluidcell.zeros(obj.N, size(e_eff,2), obj.Ntypes);
            rho_old     = fluidcell.zeros(obj.N, size(e_eff,2), obj.Ntypes);
            error_rel   = 1;
            count       = 0;

            while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
                rho = (1 - kernel*(obj.rapid_w.*rho_old))./(2*pi*(1+exp(e_eff)));
                
                % calculate error
                v1          = flatten(rho);
                v2          = flatten(rho_old);

                sumeff      = sum( v1.^2 ,1);            
                error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
                rho_old     = rho;

                count       = count+1;
            end
        end % end nested function
        
    end
    
    
    function [mu0_fit, fill_fit] = fitDensity(obj, T, density_target, mu0_guess, silent)
        % =================================================================
        % Purpose : Assuming a homogeneous system, fit the chemical
        %           potential to reproduce the desired atomic density. 
        % Input :   T       -- Temeprature of system (scalar or @(x)).
        %           density_target -- Density to fit to.
        %           mu0_guess -- Initial guess for central chemical pot.
        % Output:   mu0_fit -- Fitted central chemical potential.
        %           fill_fit -- Fitted filling function.
        % =================================================================
        if nargin < 5
            silent = false;
        end
        
%         fitfunc     = @(mu0) abs( density_target - calcDens(obj, mu0, T) );
        fitfunc     = @(mu0) abs( density_target - calcDens(obj, mu0, T) )./density_target; % relative error
        
        if ~silent
            options     = optimset('Display','iter');
            mu0_fit     = fminsearch(fitfunc, mu0_guess, options);
        else
            mu0_fit     = fminsearch(fitfunc, mu0_guess);
        end
        
        [~, fill_fit] = calcDens(obj, mu0_fit, T);
        
        function [density, fill] = calcDens(obj, mu0, T)
            mu_old      = obj.couplings{1,1};
            obj.couplings{1,1} = @(t,x) mu0;
            
            ebare   = obj.getBareEnergy(0, 0, obj.rapid_grid, obj.type_grid);

            if isa(T, 'function_handle')
                T       = T(0);
            end
            
            w       = ebare./T;
            e_eff   = obj.calcEffectiveEnergy(w, 0, 0);
            fill   = obj.calcFillingFraction(e_eff);
            
            dp      = obj.getMomentumRapidDeriv(0, 0, obj.rapid_grid, obj.type_grid);
            h0      = ones(obj.N, 1);            
            h0_dr   = obj.applyDressing(h0, fill, 0);   
            density = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*fill.*h0_dr) , 3) , 1));
            
            obj.couplings{1,1} = mu_old;
        end % end nested function
    end
    

    function [x, fill_fit] = fitThermalState(obj, fill_noneq, t, x0, options)
        % =================================================================
        % Purpose : Given a non-equilibrium state, calculate a thermal
        %           state with the same energy.
        % Input :   fill_noneq -- Filling fraction at time t.
        %           t           -- Time.
        %           x0          -- Initial guesses:
        %                           x0(1) = temperature
        %                           x0(2) = central chemical potential
        %           options     -- Fitting options.
        % Output:   x           -- Fitted temperature and mu0.
        %           fill_fit   -- Fitted therml state.
        % =================================================================
                
        if nargin < 6 % No options supplied
            options = [];
        end
        
        [N_target, E_target] = computeNE(fill_noneq, t);

        xLast   = []; % Last place computeall was called
        Na_last = []; % Use for N_atoms at xLast
        Et_last = []; % Use for E_total at xLast
        fill_last = []; 

        fun     = @objfun; % The objective function, nested below
        cfun    = @constr; % The constraint function, nested below

        % Call fmincon
        [x,f,eflag,outpt] = fmincon(fun,x0,[],[],[],[],[],[],cfun,options);
        [~,~,fill_fit] = computeall(x);

        function f = objfun(x)
            if ~isequal(x,xLast) % Check if computation is necessary
                [Na_last, Et_last, fill_last] = computeall(x);
                xLast = x;
            end
            
            % Now compute objective function
            f = abs( E_target - Et_last );
        end

        function [c,ceq] = constr(x)
            if ~isequal(x,xLast) % Check if computation is necessary
                [Na_last, Et_last, fill_last] = computeall(x);
                xLast = x;
            end
            % Now compute constraint function
            c   = [];
            ceq = N_target - Na_last;
        end
        
        function [N_atoms, E_total, fill] = computeall(x) 
            T           = x(1);
            mu0         = x(2);
            
            % Calculate thermal state
            mu0_old     = obj.couplings{1,1}(t,0);
            mu_old      = obj.couplings{1,1};
            mu_fit      = @(t,x) mu0 - (mu0_old - mu_old(t,x));
            obj.couplings{1,1} = mu_fit;
            
            ebare   = obj.getBareEnergy(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            if isa(T, 'function_handle')
                T       = T(obj.x_grid);
            end

            w       = ebare./T;
            e_eff   = obj.calcEffectiveEnergy(w, t, obj.x_grid);
            fill = obj.calcFillingFraction(e_eff);

            obj.couplings{1,1} = mu_old;
            
            [N_atoms, E_total] = computeNE(fill, t);
        end
        
        function [N_atoms, E_total] = computeNE(fill, t) 
            % Calculate number of atoms and total energy
            dp      = obj.getMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);

            h0      = obj.getOneParticleEV( 0, t, obj.x_grid, obj.rapid_grid);           
            h0_dr   = obj.applyDressing(h0, fill, t);
            
            h2      = obj.getOneParticleEV( 2, t, obj.x_grid, obj.rapid_grid);           
            h2_dr   = obj.applyDressing(h2, fill, t);
                
            ndens   = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*fill.*h0_dr) , 3) , 1));
            edens   = 1/(2*pi) * squeeze(sum( obj.rapid_w .* sum( double(dp.*fill.*h2_dr) , 3) , 1));
            
            N_atoms = trapz(obj.x_grid, ndens);
            E_total = trapz(obj.x_grid, edens);
        end  % end nested functions
          
    end
    
    
    function [v_eff, a_eff, de_dr, dp_dr] = calcVelocitiesNormal(obj, fill, t, x, rapid, type)        
        % =================================================================
        % Purpose : Overloads superclass method, as acceleration in LL
        %           model can be computed in more efficient way.
        % Input :   fill  -- filling function (fluidtensor)
        %           t     -- time (scalar)
        %           x     -- x-coordinate (can be scalar or vector)
        %           rapid -- rapid-coordinate (can be scalar or vector)
        %           type  -- quasiparticle type (can be scalar or vector)
        % Output:   v_eff -- Effective velocity (fluidtensor)
        %           a_eff -- Effective acceleration (fluidtensor)
        % =================================================================
        
        de_dr   = obj.applyDressing(obj.getEnergyRapidDeriv(t, x, rapid, type), fill, t);
        dp_dr   = obj.applyDressing(obj.getMomentumRapidDeriv(t, x, rapid, type), fill, t);
        
        v_eff   = de_dr./dp_dr;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = fluidcell.zeros( size(v_eff) );
            return
        end
        
        % Calculate acceleration from inhomogenous potential. Note dmudt
        % does not contribute as f = 0;
        a_eff_mu = 0;
        if ~isempty(obj.couplings{3,1})
            a_eff_mu = obj.couplings{3,1}(t,x);
            if size(a_eff_mu,1) == 1
                a_eff_mu = repmat(a_eff_mu, length(rapid), 1); 
            end
        end
        a_eff_mu = fluidcell(a_eff_mu);
        
        % Calculate acceleration from inhomogenous interaction
        a_eff_c = 0;
        if ~isempty(obj.couplings{2,2}) || ~isempty(obj.couplings{3,2})
            % Calculate derivative of scattering phase with respect to
            % interaction c           
            dT      = obj.getScatteringCouplingDeriv(2, t, x, rapid, obj.rapid_aux, type, obj.type_aux);
            B       = 1/(2*pi) * dT.*transpose(obj.rapid_w .* fill);
        end
        
        if ~isempty(obj.couplings{2,2}) % calc time deriv contribution
            f       = B*dp_dr;
            f_dr    = obj.applyDressing(f, fill, t);
            a_eff_c = a_eff_c + obj.couplings{2,2}(t,x).*f_dr;
        end

        if ~isempty(obj.couplings{3,2}) % calc space deriv contribution
            L       = B*de_dr;
            L_dr    = obj.applyDressing(L, fill, t);
            a_eff_c = a_eff_c + obj.couplings{3,2}(t,x).*L_dr;
        end
        
        a_eff_c = a_eff_c./dp_dr;
        a_eff   = a_eff_c + a_eff_mu;
    end
    
    
    function [v_eff, a_eff, de_dr, dp_dr] = calcVelocitiesFast(obj, fill, t, D)        
        % =================================================================
        % Purpose : Calculates effective velocity and acceleration of
        %           quasiparticles.
        % Input :   fill  -- filling function (fluidtensor)
        %           t     -- time (scalar)
        %           D     -- dressing operator (fluidtensor)
        % Output:   v_eff -- Effective velocity (fluidtensor)
        %           a_eff -- Effective acceleration (fluidtensor)
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
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = fluidcell.zeros( size(v_eff) );
            return
        end
        
        % Calculate acceleration from inhomogenous potential. Note dmudt
        % does not contribute as f = 0;
        a_eff_mu = 0;
        if ~isempty(obj.couplings{3,1})
            a_eff_mu = obj.couplings{3,1}(t,x);
            if size(a_eff_mu,1) == 1
                a_eff_mu = repmat(a_eff_mu, length(rapid), 1); 
            end
        end
        a_eff_mu = fluidcell(a_eff_mu);
        
        % Calculate acceleration from inhomogenous interaction
        a_eff_c = 0;
        if ~isempty(obj.couplings{2,2}) || ~isempty(obj.couplings{3,2})
            % Calculate derivative of scattering phase with respect to
            % interaction c           
            dT      = obj.getScatteringCouplingDeriv(2, t, x, rapid, obj.rapid_aux, type, obj.type_aux);
            B       = 1/(2*pi) * dT.*transpose(obj.rapid_w .* fill);
        end
        
        if ~isempty(obj.couplings{2,2}) % calc time deriv contribution
            f       = B*dp_dr;
            f_dr    = obj.dress(f, D);
            a_eff_c = a_eff_c + obj.couplings{2,2}(t,x).*f_dr;
        end

        if ~isempty(obj.couplings{3,2}) % calc space deriv contribution
            L       = B*de_dr;
            L_dr    = obj.dress(L, D);
            a_eff_c = a_eff_c + obj.couplings{3,2}(t,x).*L_dr;
        end
        
        a_eff_c = a_eff_c./dp_dr;
        a_eff   = a_eff_c + a_eff_mu;
    end
    
    
    function g_n = calcLocalCorrelator(obj, n, fill, t_array)
        % =================================================================
        % Purpose : Calculates local n-body correlation function
        % Input :   n     -- order of correlation calculated
        %           fill  -- filling function (cell array of iFluidTensor)
        %           t_array- array of times corresponding to fill
        % Output:   g_n   -- Correlation func. for each time in t_array
        % =================================================================
        Nsteps = length(t_array);
        g_n = zeros(obj.M, Nsteps);
        
        for k = 1:Nsteps
            if Nsteps == 1
                fill_k = fill;
                t       = t_array;
            else
                fill_k = fill{k};
                t       = t_array(k);
            end
            
            D       = obj.calcCharges(0, fill_k, t)'; % density
            prefac  = factorial(n)^2 * (obj.couplings{1,2}(t,obj.x_grid)).^n / 2^n ./ D.^n;
        
            % Find integers m_j to sum over
            m_seq   = findMseq(1, zeros(1,n), n, []);
        
            % Calc B functions
            B = obj.calcB(n, fill_k, t);
            
            g_temp = 0;
            for i = 1:length(m_seq) % for each set {m_j}_i

                prod_temp = 1;
                for j = 1:n % for each member m_j in set
                    m_j     = m_seq{i}(j);

                    prod_temp = prod_temp .* ( 1/factorial(m_j) * (B{j}./(pi*obj.couplings{1,2}(t,obj.x_grid))).^m_j );
                end

                g_temp = g_temp + prod_temp;
            end
            
            g_n(:,k) = prefac.*squeeze(double(g_temp));
        end
        
        % Nested recursive function for finding all sets of {m_j}, which 
        % satisfy the equation:
        %   sum_{j = 1}^{imfty} j*m_j = n
        % where m_j are positive integers 
        function m_seq = findMseq(j, m, n, m_seq)
            % Set first value to max possible
            R   = n - m*(1:n)';
            m(j) = ceil(R/j);

            while m(j) >= 0
                % Calculate remainder
                R_cur = n - m*(1:n)';

                if R_cur == 0
                    % save this configuration
                    m_seq{end+1} = m;
                elseif R_cur < 0
                    % no configurations possible
                else % remainder is positive
                    % look for more solutions later in sequence, but not further
                    % than the sequence length!
                    if j == n; return; end
                    m_seq = findMseq( j+1, m, n, m_seq);
                end

                m(j) = m(j) - 1;
            end
        end % end nested function
        
    end
    
        
    function B = calcB(obj, n, fill, t)
        % =================================================================
        % Purpose : Supporting function for calcLocalCorrelator()
        % Input :   n     -- order of correlation calculated
        %           fill  -- filling function (single iFluidTensor)
        %           t     -- time 
        % Output:   B     -- cell array of all B_i funcs up to order 2*n-1
        % =================================================================
        b       = cell( 1 , 2*n - 1 + 2); % added two dummy indices
        b(:)    = {fluidcell.zeros( obj.N , obj.M)};
        
        kernel1 = -1/(2*pi)*obj.getScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_aux, obj.type_grid, obj.type_aux);
        kernel2 = -(obj.rapid_grid - permute(obj.rapid_grid, [4 2 3 1])).*kernel1./obj.couplings{1,2}(t,obj.x_grid);
        
        I       = fluidcell.eye(obj.N, obj.Ntypes);

        
        X1      = I - kernel1.*transpose(obj.rapid_w.*fill);
        
        for i = 1:(2*n - 1)                
            if mod(i,2) == 0 % i even
                X2      = -kernel1*(obj.rapid_w.*fill.*b{i-2+2}) + kernel2*(obj.rapid_w.*fill.*( 2*b{i-1+2} - b{i-3+2} ));
                b{i+2}  = X1\X2;
            else % i odd
                X2      = kernel2*(obj.rapid_w.*fill.*b{i-1+2}) - kernel1*(obj.rapid_w.*fill.*b{i-2+2});
                
                if i == 1 % delta function contribution for n = 0
                    X2 = X2 + 1;
                end
                
                b{i+2}  = X1\X2;
            end
        end
        
        B = cell(1, n);
        for i = 1:n
            B{i} = 1/i*transpose(fill)*(obj.rapid_w.*b{2*i - 1 + 2 });
        end
    end
    
    
    function P = calcExcitationProb(obj, t, x, k, q)
        c = obj.couplings{1,2}(t,x);
        P = k.*q.*c.^2 ./(k.^2 .* q.^2 + 0.25*c.^2 .* (k + q).^2);
        
        P(isnan(P)) = 0; 
    end
    
    function P = calcExchangeProb(obj, t, x, rapid1, rapid2)
        c = obj.couplings{1,2}(t,x);
        P = (c.^2 .* abs(rapid1 - rapid2))./(c.^2 + (rapid1 - rapid2).^2);
        
        P(isnan(P)) = 0; 
    end
    
    
    function F = calcBackflow(obj, fill, t, x)
        phi     = -2*atan( (obj.rapid_grid - obj.rapid_aux)./obj.couplings{1,2}(t,x) );
        F       = obj.applyDressing( phi, fill, t )/2/pi;
    end
    
    
    
    
    
      
end % end public methods

    
end % end classdef