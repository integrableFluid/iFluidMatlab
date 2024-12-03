function contour_t = propagate_contour_constpiston_linear(K0, t_array, vp, veff_func, M_DSW)
    % =====================================================================
    % Purpose : Propagate Fermi contour under the action of a constant
    %           velocity piston using Matlab ODE45. 
    %           Piston moves from left to right and starts at position x=0.
    % Input :   K0        -- Fermi momentum of initial contour/background
    %           t_array   -- Vector of evolution times
    %           vp        -- Piston velocity
    %           veff_func -- Anonymous function to calculate effective
    %                           velocity for contour
    %           M_init    -- Number of points in DSW edge (R3) 
    % Output:   contour_t -- Cellarray with contour for each evolution time
    % =====================================================================


    if vp < 2*K0
        contour = [flip(linspace(0, eps, M_DSW))', flip(linspace(K0, K0+vp, M_DSW))'];
    else
        contour = [flip(linspace(0, eps, M_DSW))', flip(linspace(vp-K0, vp+K0, M_DSW))'];
    end

    N_steps     = length(t_array) - 1;
    contour_t   = cell(1, length(t_array));
    contour_t{1}= contour;

    zmin        = -eps;
    zmax        = 2*(vp+K0)*t_array(end);
    
    % initialize progress bar
    cpb = ConsoleProgressBar();                 % create instance
    initialize_progbar;                         % initialize progress bar with standard parameters
    fprintf('Time evolution progress:');
    cpb.start();   
    for i = 1:N_steps
        dt              = t_array(i+1) - t_array(i);
        contour         = step_piston(contour, dt, t_array(i));
        contour_t{i+1}  = contour;
        
        % show progress
        cpb_text = sprintf('%d/%d steps evolved', i, N_steps);
        cpb.setValue(i/N_steps);
        cpb.setText(cpb_text);
    end
    fprintf('\n')
    
    
    function contour = step_piston(contour, dt, t)
        % Gamma_pist are the points spawned at piston

        if vp < 2*K0
            cont_temp   = [ contour;     
                            zmax, K0; 
                            zmax, -K0;
                            zmin, -K0;
                            zmin, K0 + vp;
                            ];
        else
            cont_temp   = [ contour;
                            vp*t, K0;
                            zmax, K0; 
                            zmax, -K0;
                            zmin, -K0;
                            zmin, K0 + vp;
                            ];
        end

        % Calculate effective velocity and take linear step
        veff        = veff_func(cont_temp);
        contour     = contour + dt*veff(1:size(contour, 1), :);

        % Generate new points on piston edge
        if vp < 2*K0
            p_new       = [vp*(t+dt), K0 + vp];
            contour     = [p_new; contour ];
        else
            p1_new      = [vp*(t+dt), K0 + vp];
            p2_new      = [vp*(t+dt), -K0 + vp];
            contour     = [p1_new; contour; p2_new ];
        end
                 
    end
end