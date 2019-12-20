classdef XXZchainModel < iFluidCore

    % First coupling is B
    % second coupling is theta = acosh(Delta)
properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = XXZchainModel(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options)   
        % Set default parameters
        if nargin < 6
            Options = struct;
        end
   
        % Call superclass constructor
        obj = obj@iFluidCore(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
    end
    
    
    %% Implementation of abstract equations
    function ebare = getBareEnergy(obj, t, x, rapid, type)
        ebare = -sinh((obj.couplings{1,2}(t,x))).*sinh(type.*(obj.couplings{1,2}(t,x)))./(cosh(type.*(obj.couplings{1,2}(t,x)))-cos(2.*rapid))-type.*(obj.couplings{1,1}(t,x));
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        pbare = 2.*atan(coth(type.*(obj.couplings{1,2}(t,x))./2).*tan(rapid));
    end
    
    
    function de = calcEnergyRapidDeriv(obj, t, x, rapid, type)
        de = (2.*sin(2.*rapid).*sinh((obj.couplings{1,2}(t,x)).*type).*sinh((obj.couplings{1,2}(t,x))))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type)).^2;
    end

    
    function dp = calcMomentumRapidDeriv(obj, t, x, rapid, type)
        dp = -(2.*sinh((obj.couplings{1,2}(t,x)).*type))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type));
    end
    
    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Reshape input to right dimensions        
        rapid1  = reshape(rapid1, max(size(rapid1,1), size(rapid1,2)), 1, max(size(rapid1,3), size(rapid1,4)), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, max(size(rapid2,1), size(rapid2,2)), 1, max(size(rapid2,3), size(rapid2,4))); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.calcMomentumRapidDeriv(t, x, r_arg, abs(type1-type2));
        dT2     = obj.calcMomentumRapidDeriv(t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        % Calculate 3rd term
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    r_arg_temp = r_arg(:, :, min(i, size(r_arg,3)), min(j ,size(r_arg,4)) );
                    
                    temp = 2*obj.calcMomentumRapidDeriv(t, x, r_arg_temp, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = iFluidTensor( dT );
    end
    
    
    function de = calcEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. B
            de = repmat(-type, length(rapid), 1);
        else
            % Derivative w.r.t theta
            de = -(type.*sinh((obj.couplings{1,2}(t,x)))+cosh((obj.couplings{1,2}(t,x)).*type).*(sinh((obj.couplings{1,2}(t,x)).*type).*cosh((obj.couplings{1,2}(t,x))) ...
                 -type.*cos(2.*rapid).*sinh((obj.couplings{1,2}(t,x))))-cos(2.*rapid).*sinh((obj.couplings{1,2}(t,x)).*type).*cosh((obj.couplings{1,2}(t,x))))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type)).^2;
        end
    end

    
    function dp = calcMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. B
            dp = 0;
        else
            % Derivative w.r.t theta
            dp = (type.*sin(2.*rapid))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type));
        end
    end
    
    
    function dT = calcScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        if coupIdx == 1 % deriv w.r.t. B
            dT  = 0;
            dT  = iFluidTensor( dT );
            return
        end 
        % Else it's deriv w.r.t. theta
        
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.calcMomentumCouplingDeriv(2, t, x, r_arg, abs(type1-type2));
        dT2     = obj.calcMomentumCouplingDeriv(2, t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        % Calculate 3rd term
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.calcMomentumCouplingDeriv(2, t, x, r_arg, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = iFluidTensor( dT );
    end
    
    
    
    
      
end % end public methods

    
end % end classdef