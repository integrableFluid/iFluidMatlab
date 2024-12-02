function [G_t, G_top_t, G_bot_t] = propagate_contour_periodic(G_top, G_bot, t_array, L, veff_func)
    
    % Time evolution of contour with periodic boundary conditions confined
    % to region [-L/2, L/2]. Uses linear stepping.
    % To enforce PBC, it is easier working with a disconnected contour, 
    % split into a top and bottom part. These parts are generally
    % undordered, meaning that the first index does not correspond to the
    % first point in the contour.


    % Check if the left and right side points overlap once periodic BC are
    % enforced. If so, delete the end-point. 
    % Otherwise the find_crossings function will cause an error!
    if G_top(1, 1) == -L/2 && G_top(end, 1) == L/2
        G_top(end, :) = []; 
    end
    if G_bot(1, 1) == L/2 && G_bot(end, 1) == -L/2
        G_bot(end, :) = []; 
    end

    Nsteps      = length(t_array) - 1;
    G_top_t     = cell(1, length(t_array));
    G_bot_t     = cell(1, length(t_array));
    G_t         = cell(1, length(t_array));

    G_top_t{1}  = G_top; 
    G_bot_t{1}  = G_bot;

    % initialize progress bar
    cpb = ConsoleProgressBar();                 % create instance
    initialize_progbar;                         % initialize progress bar with standard parameters
    fprintf('Time evolution progress:');
    cpb.start();   
    for i = 1:Nsteps

        % combine top and bottom contours and add periodic padding
        [G_tot, it, ib] = merge_contours_periodic(G_top, G_bot, L);
        G_t{i}          = G_tot;

        % calculate effective velocity
        veff            = veff_func(G_tot);

        % take linear step
        dt              = t_array(i+1) - t_array(i);
        G_tot           = G_tot + dt*veff;

        % separate contour back into top and bottom part
        G_top           = G_tot(it, :);
        G_bot           = G_tot(ib, :);
        
        % enforce periodic boundary conditions       
        G_bot(G_bot(:,1) < -L/2, 1)  = G_bot(G_bot(:,1) < -L/2, 1) + L;
        G_top(G_top(:,1) > L/2, 1)  = G_top(G_top(:,1) > L/2, 1) - L;

        % store the contours
        G_top_t{i+1}    = G_top;
        G_bot_t{i+1}    = G_bot;
        
        % show progress
        cpb_text = sprintf('%d/%d steps evolved', i, Nsteps);
        cpb.setValue(i/Nsteps);
        cpb.setText(cpb_text);
    end
    fprintf('\n')

    % store merged contour at final time step
    G_tot       = merge_contours_periodic(G_top, G_bot, L);
    G_t{i+1}    = G_tot;
    
end