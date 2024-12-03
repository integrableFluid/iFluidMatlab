function Gamma_t = propagate_contour_linear(Gamma, t_array, veff_func)
    % =====================================================================
    % Purpose : Propagate Fermi contour with linear step.
    % Input :   Gamma       -- initial contour with each row = (z, rapid)
    %           t_array     -- array of evolution times
    %           veff_func   -- anonymous function to calculate effective
    %                           velocity (and acceleration) for contour
    % Output:   Gamma_t     -- cell array with contour for each time
    % =====================================================================

    Nsteps      = length(t_array) - 1;
    Gamma_t     = cell(1, length(t_array));
    Gamma_t{1}  = Gamma; 

  
    % initialize progress bar
    cpb = ConsoleProgressBar();                 % create instance
    initialize_progbar;                         % initialize progress bar with standard parameters
    fprintf('Time evolution progress:');
    cpb.start();   
    for i = 1:Nsteps

        % calculate effective velocity
        veff        = veff_func(Gamma);

        % take linear step
        dt          = t_array(i+1) - t_array(i);
        Gamma       = Gamma + dt*veff;

        % store the contours
        Gamma_t{i+1}= Gamma;
        
        % show progress
        cpb_text = sprintf('%d/%d steps evolved', i, Nsteps);
        cpb.setValue(i/Nsteps);
        cpb.setText(cpb_text);
    end
    fprintf('\n')
end