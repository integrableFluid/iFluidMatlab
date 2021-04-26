classdef XXZchainModel < iFluidCore
    % iFluid implmentation of TBA functions for XXZ chain model
    %
    % ### Couplings are { B, acosh(Delta) } ###
    %
    % NOTE: This TBA is only valid for B < 0 and Delta > 1
    %
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
    
    
    function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
        de = (2.*sin(2.*rapid).*sinh((obj.couplings{1,2}(t,x)).*type).*sinh((obj.couplings{1,2}(t,x))))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type)).^2;
    end

    
    function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
        dp = -(2.*sinh((obj.couplings{1,2}(t,x)).*type))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type));
    end
    
    
    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1, 1);
        I_type  = permute(I_type, [3 5 1 4 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.getMomentumRapidDeriv(t, x, r_arg, abs(type1-type2));
        dT2     = obj.getMomentumRapidDeriv(t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        % Calculate 3rd term
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    r_arg_temp = r_arg(:, :, min(i, size(r_arg,3)), min(j ,size(r_arg,5)) );
                    
                    temp = 2*obj.getMomentumRapidDeriv(t, x, r_arg_temp, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,:,j) = dT3(:,:,i,:,j) + temp;
                end
            end
        end
        
        dT  = fluidcell( dT1 + dT2 + dT3 );
    end
    
    
    function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. B
            de = repmat(-type, length(rapid), 1);
        else
            % Derivative w.r.t theta
            de = -(type.*sinh((obj.couplings{1,2}(t,x)))+cosh((obj.couplings{1,2}(t,x)).*type).*(sinh((obj.couplings{1,2}(t,x)).*type).*cosh((obj.couplings{1,2}(t,x))) ...
                 -type.*cos(2.*rapid).*sinh((obj.couplings{1,2}(t,x))))-cos(2.*rapid).*sinh((obj.couplings{1,2}(t,x)).*type).*cosh((obj.couplings{1,2}(t,x))))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type)).^2;
        end
    end

    
    function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 1
            % Derivative w.r.t. B
            dp = 0;
        else
            % Derivative w.r.t theta
            dp = (type.*sin(2.*rapid))./(cos(2.*rapid)-cosh((obj.couplings{1,2}(t,x)).*type));
        end
    end
    
    
    function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        if coupIdx == 1 % deriv w.r.t. B
            dT  = 0;
            dT  = fluidcell( dT );
            return
        end 
        % Else it's deriv w.r.t. theta
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1, 1);
        I_type  = permute(I_type, [3 5 1 4 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.getMomentumCouplingDeriv(2, t, x, r_arg, abs(type1-type2));
        dT2     = obj.getMomentumCouplingDeriv(2, t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        % Calculate 3rd term
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.getMomentumCouplingDeriv(2, t, x, r_arg, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,:,j) = dT3(:,:,i,:,j) + temp;
                end
            end
        end
 
        dT  = fluidcell( dT1 + dT2 + dT3 );
    end
    
    
    function h_i = getOneParticleEV(obj, charIdx, t, x, rapid)
        % Overload the iFluidCore function, because h_0 counters higher
        % types as multiple particles.
        switch charIdx
        case 0 % eigenvalue of number operator
            h_i = repmat(obj.type_grid, length(rapid), 1);
        case 1 % eigenvalue of momentum operator
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of Hamiltonian
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid);
        otherwise 
            % Higher order charges must be implemented in specific model
            error(['Eigenvalue ' num2str(charIdx) ' not implmented!'])
        end
    end
    
      
end % end public methods

    
end % end classdef