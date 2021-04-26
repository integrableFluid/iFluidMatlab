classdef HardRodModel < iFluidCore

    
properties (Access = protected)

    % Species of quasiparticle
    quasiSpecies= 'classical'; 
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = HardRodModel(x_grid, rapid_grid, rapid_w, couplings, Options)   
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
        ebare = rapid.^2 /2 - obj.couplings{1,2}(t,x);
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        pbare = rapid;
    end
    
    
    function de = getEnergyRapidDeriv(obj, t, x, rapid, type)
        de = rapid;
    end

    
    function dp = getMomentumRapidDeriv(obj, t, x, rapid, type)
        dp = repmat(1, length(rapid), 1);
    end
    
    
    function dT = getScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        dT      = repmat(obj.couplings{1,1}(t,x) , length(rapid1), 1 , 1 , length(rapid2));         
        dT      = fluidcell(dT); % Converts to iFluidTensor
    end
    
    
    function de = getEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        if coupIdx == 2
            de = repmat(-1, length(rapid), 1);
        else
            de = 0;
        end
    end

    
    function dp = getMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
       dp = 0;
    end
    
    
    function dT = getScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        if coupIdx == 1
            dT = repmat(1 , length(rapid1), 1 , 1 , length(rapid2));
        else
            dT = 0;
        end

        dT = fluidcell(dT); % Converts to iFluidTensor
    end
    
      
end % end public methods

    
end % end classdef