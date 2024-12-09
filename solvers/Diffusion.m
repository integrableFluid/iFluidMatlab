classdef Diffusion < AdvectionSource
    
    
properties (Access = protected)
    
    
end % end protected properties
    

methods (Access = public)
    
    % Superclass constructor
    function obj = Diffusion(model, varargin)        
        obj = obj@AdvectionSource(model, varargin{:});
        

    end
    

    function [D, DT] = calcDiffusion(obj, fill, t)
        % =================================================================
        % Purpose : Calculate diffusion kernel from given filling fill.
        % Input :   fill   -- Filling function at time t.
        %           t       -- Time.
        % Output:   D       -- Diffusion matrix.
        %           DT      -- Diffusion kernel.
        % =================================================================
                
        rapid_aux   = permute(obj.rapid_grid, [4 2 3 1]);
        type_aux    = permute(obj.type_grid, [1 2 5 4 3]);
        [rhoP, rhoS]= obj.model.transform2rho(fill, t);
        
        % Prepare required quantities
        I           = fluidcell.eye(obj.N, obj.Ntypes);
        T           = -1/(2*pi)*obj.model.getScatteringRapidDeriv(0, obj.x_grid, obj.rapid_grid, rapid_aux, obj.type_grid, type_aux);
        T_dr        = obj.model.applyDressing(T, fill, 0);
        v_eff       = obj.model.calcEffectiveVelocities(fill, 0);

            
        % Calculate diffusion kernel
        W           = rhoP.*(1-fill).*T_dr.^2 .* abs(v_eff - v_eff.t());
        w           = sum( W.*obj.rapid_w, 1); 
        DT          = rhoS.^(-2).*(I.*w./obj.rapid_w - W).*sqrt(obj.rapid_w.*permute(obj.rapid_w, [4 2 3 1]));
        
        % Calculate diffusion operator      
        R           = (I - T.*(fill.*obj.rapid_w))./rhoS;   
        fill_x      = obj.gradient_space(obj.x_grid, fill);
        temp        = double( R\(DT*fill_x) );   
        D           = 0.5*R*obj.gradient_space(obj.x_grid, temp);
        
    end


    function f_x = gradient_space(obj, x, f)
        % Calculate spatial gradient of tensor f
        
        f       = double(f);
        dx      = x(2) - x(1);
        
        % repeat endpoints for circular deriv
        f       = [f(:,1,:), f, f(:,end,:)]; 
        
        % dnumerical derivative via central difference 
        f_x     = (circshift(f,-1,2)-circshift(f,1,2))/(2*dx);
        
        % remove padding
        f_x     = f_x(:,2:end-1);

        % ensure periodicity
        if strcmpi(obj.SolverObj.getSettings.boundary_conditions, 'periodic')
            f_edge      = 0.5*(f_x(:,1,:) + f_x(:,end,:));
            
            f_x(:,1,:)  = 2*f_edge;
            f_x(:,end,:)= 2*f_edge;
        end
        
        f_x         = fluidcell(f_x);
    end
    
    
    % Implementation of abstract functions
    
    function [fill, u, w] = initialize(obj, fill_init, u_init, w_init, t_array)
        % =================================================================
        % Purpose : Calculates and stores all quantities needed for the
        %           step() method.
        %           In this case of first order step, nothing is required.
        % Input :   fill_init -- Initial filling function.
        %           u_init     -- Initial auxiliary variable 1.
        %           w_init     -- Initial auxiliary variable 2.
        %           t_array    -- Array of time steps.
        % Output:   fill      -- Input fill for first step().
        %           u          -- Input u for first step().
        %           w          -- Input w for first step().
        % =================================================================
       
        if iscell(fill_init) && length(fill_init) > 1
            % Assume fill_init{end} = fill(t = 0),
            % while fill_init{end-1} = fill(t = -dt), ...
            dt          = t_array(2) - t_array(1);
            
            [S, A]      = obj.calcSourceTerm(fill_init{end-1}, u_init, w_init, -dt);
            
            obj.S_prev  = S;
            obj.A_prev  = A;
            
            fill       = fill_init{end};
            u           = u_init; 
            w           = w_init; 
        else
            
            [S, A]      = obj.calcSourceTerm(fill_init, u_init, w_init, 0);
            
            obj.S_prev  = S;
            obj.A_prev  = A;
            
            fill       = fill_init;
            u           = u_init; 
            w           = w_init;
        end

    end
    
    
    function [S, A] = calcSourceTerm(obj, fill, u, w, t)
        S   = obj.calcDiffusion(fill, t);
        A   = 0; % no auxiliary source term for diffusion
    end
      
    
    function [u_next, w_next] = auxStep(obj, u, w, A, t, dt)
        % no auxiliary source term --> do nothing
        u_next = u;
        w_next = w;
    end
    
    
     
end % end




end % end classdef