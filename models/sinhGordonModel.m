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
    
    
    function de = calcEnergyRapidDeriv(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        de = m.*sinh(rapid);
    end

    
    function dp = calcMomentumRapidDeriv(obj, t, x, rapid, type)
        m = sqrt( obj.couplings{1,2}(t,x).^2 .* sin(pi*obj.couplings{1,1}(t,x))./(pi*obj.couplings{1,1}(t,x)) );
        dp = m.*cosh(rapid);
    end
    
    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Reshape input to ensure right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, 1, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, 1, length(type2)); % type2 is 4th index
        
        dT = -(2.*sin((obj.couplings{1,1}(t,x).*pi)).*cosh((rapid1-rapid2)))./(sinh((rapid1-rapid2)).^2+sin(obj.couplings{1,1}(t,x).*pi).^2);
        
        dT(isnan(dT)) = 0; % removes any NaN
        dT      = iFluidTensor(dT); % Converts to iFluidTensor
    end
    
    
    function de = calcEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
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

    
    function dp = calcMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
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
    
    
    function dT = calcScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)        
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, 1, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, 1, length(type2)); % type2 is 4th index
        
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
        dT  = iFluidTensor( dT );
    end
    
    
    function Psi_k_t = calcVertexExpval( obj, kmax, theta_t, t_array)
        % =================================================================
        % Purpose : Calculates expectation value of local vertex operator
        %           for k = 1 to kmax
        % Input :   kamx    -- max order of expval calculated
        %           theta_t -- filling function (cell array of iFluidTensor)
        %           t_array -- array of times corresponding to theta
        % Output:   Psi_k_t -- Cell array of expectation values for each
        %                      time sgiven in t_array 
        % =================================================================
        Psi_k_t = cell(1, length(t_array));

        for i = 1:length(t_array)
            t       = t_array(i);
            theta   = theta_t{i};
            Psi_mat = zeros(obj.M, kmax);
            
            Psi_k   = 1;

            for k = 1:kmax
                Xi      = obj.calcXi( k, t, theta );
                Psi_k   = Psi_k .* ( 1 + 2*sin(pi*obj.couplings{1,1}(t,obj.x_grid)*(2*k + 1))/pi .* sum( double(obj.rapid_w.*Xi.*theta.*exp(obj.rapid_grid)) , 1 ) );
                
                Psi_mat(:,k) = squeeze(Psi_k);
            end

            Psi_k_t{i} = Psi_mat;
        end
    end
      
end % end public methods

    
methods (Access = private)    

    function Xi = calcXi(obj, k_idx, t, theta)
        % Supporting function for calcVertexExpval()
        rapid_arg   = obj.rapid_grid - permute(obj.rapid_grid, [4 2 3 1]);
        chi         = obj.calcChi(k_idx, t, obj.x_grid, rapid_arg);
        
        X           = permute(eye(obj.N), [1 3 4 2]) - chi.*transpose(obj.rapid_w.*theta); 
        epsi        = iFluidTensor(exp(-obj.rapid_grid));

        Xi          = X\epsi;   
    end
    

    function chi = calcChi(obj, k, t, x, rapid)
        % Supporting function for calcVertexExpval()
        chi = 1i/(2*pi) * ( exp(-1i*2*k*obj.couplings{1,1}(t,x))./sinh( rapid + 1i*pi*obj.couplings{1,1}(t,x) ) ...
                            - exp(1i*2*k*obj.couplings{1,1}(t,x))./sinh( rapid - 1i*pi*obj.couplings{1,1}(t,x) ) );
    end

end % end private methods


end % end classdef