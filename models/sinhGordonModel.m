classdef sinhGordonModel < iFluidCore
    % iFluid implmentation of TBA functions for relativistic sinh-Gordom model
    %
    % ### Couplings are { alpha, beta, mu } ###

properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = sinhGordonModel(x_grid, rapid_grid, rapid_w, couplings, Options)   
        % Set default parameters
        if nargin < 5
            Options = struct;
        end
        
        % sinh-Gordon model has 1 species of quasi-particles
        Ntypes = 1;
   
        % Call superclass constructor
        obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
    end
    
    
    %% Implementation of abstract equations
    function ebare = getBareEnergy(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        ebare = m.*cosh(rapid)  - obj.couplings{1,3}(t,x);
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        pbare = m.*sinh(rapid);
    end
    
    
    function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        de = m.*sinh(rapid);
    end

    
    function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        dp = m.*cosh(rapid);
    end
    
    
    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        dT = -(2.*sin((obj.couplings{1,1}(t,x).*pi)).*cosh((rapid1-rapid2)))./(sinh((rapid1-rapid2)).^2+sin(obj.couplings{1,1}(t,x).*pi).^2);
        
        dT(isnan(dT)) = 0; % removes any NaN
        dT      = fluidcell(dT); % Converts to iFluidTensor
    end
    
    
    function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. alpha
            de = ((pi*obj.couplings{1,1}(t,x).*cot(pi*obj.couplings{1,1}(t,x)) - 1).*sqrt((obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x)))./obj.couplings{1,1}(t,x)))./(2*sqrt(pi)*obj.couplings{1,1}(t,x)).*cosh(rapid);
        elseif coupIdx == 2
            % Derivative w.r.t mu
            de = sqrt((obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x)))./obj.couplings{1,1}(t,x))./(sqrt(pi*obj.couplings{1,2}(t,x))).*cosh(rapid);
        else
            % Derivative w.r.t chamical potential
            de = repmat(-1, length(rapid), 1);
        end
    end

    
    function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. alpha
            dp = ((pi*obj.couplings{1,1}(t,x).*cot(pi*obj.couplings{1,1}(t,x)) - 1).*sqrt((obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x)))./obj.couplings{1,1}(t,x)))./(2*sqrt(pi)*obj.couplings{1,1}(t,x)).*sinh(rapid);
        elseif coupIdx == 2
            % Derivative w.r.t mu
            dp = sqrt((obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x)))./obj.couplings{1,1}(t,x))./(sqrt(pi*obj.couplings{1,2}(t,x))).*sinh(rapid);
        else
            % Derivative w.r.t chamical potential
            dp = 0;
        end
    end
    
    
    function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)                
        if coupIdx == 1
            % Derivative w.r.t. alpha
            dT = (pi.*cos(((obj.couplings{1,1}(t,x)).*pi)).*sinh((rapid1-rapid2)))./(cosh((rapid1-rapid2)).^2-cos(((obj.couplings{1,1}(t,x)).*pi)).^2);
        elseif coupIdx == 2
            % Derivative w.r.t mu
            dT = 0;
        else
            % Derivative w.r.t chemical potential
            dT = 0;
        end
        
        dT(isnan(dT)) = 0; % removes any NaN
        dT  = fluidcell( dT );
    end
    
    
    function [Psi, V_Psi] = calcVertexExpval( obj, kmax, theta_t, t_array, calcV)
        % =================================================================
        % Purpose : Calculates expectation value of local vertex operator
        %           for k = 1 to kmax
        % Input :   kamx    -- max order of expval calculated
        %           theta_t -- filling function (cell array of iFluidTensor)
        %           t_array -- array of times corresponding to theta
        %           calcV   -- if true, calc V
        % Output:   Psi     -- vertex expactation values
        %           V_Psi   -- corresponding V-field 
        % =================================================================
        
        if nargin < 5
            calcV = false;
        end
        
        Psi = zeros( obj.M, kmax, length(t_array) );
        V_Psi = cell(kmax, length(t_array));
        
        if calcV
            [~, rhoS_t] = obj.transform2rho(theta_t, t_array);
        end

        for i = 1:length(t_array)
            t       = t_array(i);
            theta   = theta_t{i};

            Hk      = zeros( obj.M, kmax ); 
            Xi_m    = zeros( obj.N, obj.M, kmax );
            Xi_p    = zeros( obj.N, obj.M, kmax );

            for k = 1:kmax
                Hk(:,k) = obj.calcH(k-1, t, theta);
                
                if calcV
                    Xi_m(:,:,k) = double( obj.calcEps(k-1, t, theta, -1) );
                    Xi_p(:,:,k) = double( obj.calcEps(k-1, t, theta, +1) );
                end
            end
            
            Psi(:,:,i) = cumprod(Hk, 2);
            
            if calcV
                prefac = 2*sin(pi*obj.couplings{1,1}(t,obj.x_grid)*(2*(0:kmax-1) + 1))/pi;
                V_temp = cumsum( permute(prefac./Hk, [3 1 2]).*Xi_m.*Xi_p , 3 ).*permute(Psi(:,:,i),[3 1 2])./double(rhoS_t{i});
                for k = 1:kmax
                    V_Psi{k,i} = fluidcell( V_temp(:,:,k) );
                end
                
            end
        end
    end
      
end % end public methods

    
methods (Access = private)
    
    function Hk = calcH(obj, k, t, theta )
        % Supporting function for calcVertexExpval()
        eps     = obj.calcEps( k, t, theta, -1);
        Hk      = ( 1 + 2*sin(pi*obj.couplings{1,1}(t,obj.x_grid)*(2*k + 1))/pi .* sum( double(obj.rapid_w.*eps.*theta.*exp(obj.rapid_grid)) , 1 ) );
    end

    function eps = calcEps(obj, k, t, theta, sgn)
        % Supporting function for calcVertexExpval()
        rapid_arg   = - (obj.rapid_grid - permute(obj.rapid_grid, [4 2 3 1]));        
        chi         = 1/pi*imag(exp(2*k*1i*pi*obj.couplings{1,1}(t,obj.x_grid))./sinh(sign(sgn)*rapid_arg - 1i*pi*obj.couplings{1,1}(t,obj.x_grid)));
        
        X           = permute(eye(obj.N), [1 3 4 2]) - chi.*transpose(obj.rapid_w.*theta); 
        epsi        = fluidcell(exp(sign(sgn)*obj.rapid_grid));

        eps         = X\epsi;   
    end
    

end % end private methods


end % end classdef