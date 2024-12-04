function contour_t = propagate_contour_piston_linear(contour, t_array, K0, xp_func, vp_func, veff_func)
    % =====================================================================
    % Purpose : Propagate Fermi contour under the action of a constant
    %           velocity piston using linear step updates.
    %           The output contour contains only points that have been in
    %           contact with the piston and points in the initial input
    %           contour.
    %           Otherwise, it is assumed that the piston moves though a
    %           homogeneous background Fermi sea [-K0; K0], and that the 
    %           piston moves from left to right.
    % Input :   contour   -- Initial contour with each row = (z, rapid)
    %           t_array   -- Vector of evolution times
    %           K0        -- Fermi momentum/rapidity of background state
    %           xp_func   -- Piston position as anonynous function of time
    %           vp_func   -- Piston velocity as anonynous function of time
    %           veff_func -- Anonymous function to calculate effective
    %                           velocity for contour
    % Output:   contour_t -- Cellarray with contour for each evolution time
    % =====================================================================

    if isempty(contour)
        contour = [xp_func(0), K0];
    end

    N_steps     = length(t_array) - 1;
    contour_t   = cell(1, length(t_array));
    contour_t{1}= contour;

    % Calculate speed of sound for background
    vsound = veff_func([0 K0; 10 K0; 10 -K0; 0 -K0]);
    vsound = vsound(1);

    
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

        if vp_func(t) < 2*vsound
            cont_temp   = [ contour;     
                            1e10, K0; 
                            1e10, -K0;
                            -1e10, -K0;
                            -1e10, K0 + vp_func(t);
                            ];
        else
            cont_temp   = [ contour;
                            xp_func(t), K0;
                            1e10, K0; 
                            1e10, -K0;
                            -1e10, -K0;
                            -1e10, K0 + vp_func(t);
                            ];
        end

        % Calculate effective velocity and take linear step
        veff        = veff_func(cont_temp);
        contour     = contour + dt*veff(1:size(contour, 1), :);

        % Generate new points on piston edge
        if vp_func(t) < 2*vsound
            p_new       = [xp_func(t+dt), K0 + vp_func(t+dt)];
            contour     = [p_new; contour ];
        else
            p1_new      = [xp_func(t+dt), K0 + vp_func(t+dt)];
            p2_new      = [xp_func(t+dt), -K0 + vp_func(t+dt)];
            contour     = [p1_new; contour; p2_new ];
        end
                 
    end
end