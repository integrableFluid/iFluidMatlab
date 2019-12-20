classdef LiebLinigerModel_SI
    % Wrapper class for LiebLinigerModel.
    % Converts inputs from SI-units to iFluid internal units
    % Converts outputs to SI-units
    %
    % Scales all units to those of a quantum harmonic oscilator with
    % frequency omega_scale
    %
    
properties (Access = public)
    % Physical constants for Rb-87
    m_si        = 87*1.6605402e-27; % Rb-87 mass
    hbar_si     = 1.054571726e-34;
    kB_si       = 1.38065e-23;
    as_si       = 5.2e-9; % Rb-87 scattering length
    
    % TBA unit scales     
    Eg_si       = []; % energy scale
    Lg_si       = []; % length scale
    t_si        = []; % time scale
    T_si        = []; % temperature scale
    P_si        = []; % momentum scale
    
    LLM         = []; % LiebLinigerModel

end % end private properties


methods (Access = public)
    
    %% Constructor
    function obj = LiebLinigerModel_SI(omega_scale, x_grid, rapid_grid, rapid_w, couplings, Options)
        if nargin < 5
            % If no options, pass empty struct to LiebLiniger constructor
            Options = struct;
        end
        
        % Calculate unit scales
        obj.Eg_si   = 0.5*obj.hbar_si*omega_scale; 
        obj.Lg_si   = (sqrt(obj.m_si*omega_scale/obj.hbar_si))^(-1); 
        obj.t_si    = obj.hbar_si/obj.Eg_si; 
        obj.T_si    = obj.Eg_si/obj.kB_si;
        obj.P_si    = obj.hbar_si/obj.Lg_si;
        
        % Convert grids and couplings to TBA units
        x_grid      = obj.convert2TBA(x_grid, 'length');
        rapid_grid  = obj.convert2TBA(rapid_grid, 'rapidity');
        rapid_w     = obj.convert2TBA(rapid_w, 'rapidity');
        couplings   = obj.convert2TBA(couplings, 'couplings');
        
        % Instantiate LiebLinigerModel
        obj.LLM     = LiebLinigerModel(x_grid, rapid_grid, rapid_w, couplings, Options);
    end
    
    
    %% Unit-convertion functions
    function quantity_tba = convert2TBA(obj, quantity_si, unit)
        switch unit
            case 'energy'
                quantity_tba = quantity_si/obj.Eg_si;
            case 'rapidity'
                quantity_tba = quantity_si*obj.Lg_si;
            case 'momentum'
                quantity_tba = quantity_si/obj.P_si;
            case 'time'
                quantity_tba = quantity_si/obj.t_si;
            case 'length'
                quantity_tba = quantity_si/obj.Lg_si;
            case 'temperature'
                quantity_tba = quantity_si/obj.T_si;
            case 'couplings'
                % Anonymous functions are in SI units and take SI
                % arguments. Thus, convert arguments to SI and output to
                % TBA units.
                quantity_tba{1,1} = @(t, x) quantity_si{1,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si; % mu is in units of energy
                quantity_tba{1,2} = @(t, x) quantity_si{1,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si; % c is in units of inverse length
                
                % time derivatives
                quantity_tba{2,1} = @(t, x) quantity_si{2,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.t_si; 
                quantity_tba{2,2} = @(t, x) quantity_si{2,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.t_si; 
                
                % space derivatives
                quantity_tba{3,1} = @(t, x) quantity_si{3,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.Lg_si; 
                quantity_tba{3,2} = @(t, x) quantity_si{3,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.Lg_si; 
                
                % Make sure empty coupling derivatives stay empty!  
                for i = 1:size(quantity_si,1)
                    for j = 1:size(quantity_si,2)
                        if isempty( quantity_si{i,j} )
                            quantity_tba{i,j} = [];
                        end
                    end
                end
            otherwise
                error('Conversion of specified unit is not implemented!')
        end 
    end
    
    
    function quantity_si = convert2SI(obj, quantity_tba, unit)
        switch unit
            case 'energy'
                quantity_si = quantity_tba*obj.Eg_si;
            case 'momentum'
                quantity_si = quantity_tba*obj.P_si;
            case 'rapidity'
                quantity_si = quantity_tba/obj.Lg_si;
            case 'time'
                quantity_si = quantity_tba*obj.t_si;
            case 'length'
                quantity_si = quantity_tba*obj.Lg_si;
            case 'temperature'
                quantity_si = quantity_tba*obj.T_si;
            case 'couplings'
                % Couplings
                quantity_si{1,1} = @(t, x) quantity_tba{1,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si; % mu is in units of energy
                quantity_si{1,2} = @(t, x) quantity_tba{1,2}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si; % c is in units of inverse length
                
                % time derivatives
                quantity_si{2,1} = @(t, x) quantity_tba{2,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si/obj.t_si; 
                quantity_si{2,2} = @(t, x) quantity_tba{2,2}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si/obj.t_si; 
                
                % space derivatives
                quantity_si{3,1} = @(t, x) quantity_tba{3,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si/obj.Lg_si; 
                quantity_si{3,2} = @(t, x) quantity_tba{3,2}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si/obj.Lg_si; 
                
                % Make sure empty coupling derivatives stay empty!  
                for i = 1:size(quantity_si,1)
                    for j = 1:size(quantity_si,2)
                        if isempty( quantity_tba{i,j} )
                            quantity_si{i,j} = [];
                        end
                    end
                end
            otherwise
                error('Conversion of specified unit is not implemented!')
        end 
    end
    
    %% Wrapper functions    
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, mu0_guess, setCouplingFlag)        
        % Convert SI --> TBA
        T       = obj.convert2TBA(T, 'temperature');
        mu0_guess = obj.convert2TBA(mu0_guess, 'energy');
        if ~isempty(V_ext)
            V_ext   = @(t, x) V_ext( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si;
        end
        
        % Run LLS function
        mu0_fit = obj.LLS.fitAtomnumber(T, V_ext, Natoms, mu0_guess, setCouplingFlag);
        
        % Convert TBA --> SI
        mu0_fit = obj.convert2SI(mu0_fit, 'energy');
    end
    
    
    function setCouplings(obj, couplings)
        % Convert SI --> TBA
        couplings = obj.convert2TBA(couplings, 'couplings');
        
        obj.LLS.setCouplings(couplings);
    end
    
    
    function couplings = getCouplings(obj)
        couplings = obj.LLS.getCouplings();
        
        % Convert TBA --> SI
        couplings = obj.convert2SI(couplings, 'couplings');
    end
    
    
    function [theta_t, u_t] = propagateTheta(obj, theta_init, t_array)
        % Convert SI --> TBA
        t_array = obj.convert2TBA(t_array, 'time');
        
        % Run LLS function
        [theta_t, u_t] = obj.LLS.propagateTheta(theta_init, t_array);
        
        % Convert TBA --> SI
%         u_t = obj.convert2SI(u_t, 'length');
    end
    
    
    function [rho, rhoS] = transform2rho(obj, theta, t_array)
        if nargin == 3       
            % Convert SI --> TBA
            t_array = obj.convert2TBA(t_array, 'time');

            % Run LLS function
            [rho, rhoS] = obj.LLS.transform2rho(theta, t_array);
        else
            % Run LLS function
            [rho, rhoS] = obj.LLS.transform2rho(theta);
        end
    end
    
    
    function [theta, rhoS] = transform2theta(obj, rho, t_array)
        % Convert SI --> TBA (for inverse quantities, use SI rather than TBA)
        rho     = obj.convert2SI( rho, 'length'); % convert 'per length' 
        rho     = obj.convert2SI( rho, 'rapidity'); % convert 'per rapidity' 
        
        if nargin == 3       
            % Convert SI --> TBA
            t_array = obj.convert2TBA(t_array, 'time');

            % Run LLS function
            [theta, rhoS] = obj.LLS.transform2theta(rho, t_array);
        else
            % Run LLS function
            [theta, rhoS] = obj.LLS.transform2theta(rho);
        end
    end
    
    
    function [q, j] = calcCharges(obj, theta, c_idx, t_array)
        % Convert SI --> TBA
        t_array = obj.convert2TBA(t_array, 'time');
        
        % Run LLS function
        [q, j] = obj.LLS.calcCharges(theta, c_idx, t_array);
        
        % Convert TBA --> SI
        % NOTE: Doesn't convert currents
        for i = length(c_idx)
            switch c_idx(i)
            case 0 % atomic density
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            case 1 % momentum density
                q(:,:,i) = obj.convert2SI(q(:,:,i), 'momentum'); 
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            case 2 %energy density
                q(:,:,i) = obj.convert2SI(q(:,:,i), 'energy'); 
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            otherwise
                disp(['Warning: No known unit of charge nr. ' num2str(c_idx(i))])
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % Convert SI --> TBA
        T = obj.convert2TBA(T, 'temperature');
        
        if nargin == 3
            TBA_couplings   = obj.convert2TBA(TBA_couplings, 'couplings');
            [theta, e_eff]  = obj.LLS.calcThermalState(T, TBA_couplings);
        else
            [theta, e_eff]  = obj.LLS.calcThermalState(T);
        end
    end
    
    
end % end public methods

    
end % end classdef