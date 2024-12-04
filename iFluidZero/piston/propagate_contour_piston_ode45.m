function [contour_t, tout, teout] = propagate_contour_piston_ode45(cont_init, tspan, K0, M_spawn, xp_func, vp_func, veff_func)
    % =====================================================================
    % Purpose : Propagate Fermi contour under the action of an piston using
    %           Matlab's built-in ODE45 solver.
    %           The output contour contains only points that have been in
    %           contact with the piston and points in the initial input
    %           contour.
    %           Otherwise, it is assumed that the piston moves though a
    %           homogeneous background Fermi sea [-K0; K0], and that the 
    %           piston moves from left to right.
    % Input :   cont_init -- Initial contour with each row = (z, rapid)
    %           tspan     -- Matlab ODE tspan
    %           K0        -- Fermi momentum/rapidity of background state
    %           M_spawn   -- Number of points at piston spawned during evol
    %           xp_func   -- Piston position as anonynous function of time
    %           vp_func   -- Piston velocity as anonynous function of time
    %           veff_func -- Anonymous function to calculate effective
    %                           velocity for contour
    % Output:   contour_t -- Cell array with contour for each time
    %           tout      -- Corresponding evolution times
    %           teout     -- Array of times of triggered events
    % =====================================================================
      

    % For ODE solver, the number of contour points must be constant. Hence,
    % when spawning new points at piston, a leading "dummy" point must be
    % deleted; the starting contour thus contains M additional NaN points.
    if isempty(cont_init)
        cont_init = [xp_func(0), K0];
    end

    contour = [cont_init; nan(M_spawn, 2)];
    y0      = contour(:);

    Nc      = size(cont_init, 3); % number of contours
    Np      = size(cont_init, 1) + M_spawn; % number of discretization points
   
    
    options = odeset('Events',@events);

    i_spawn = 1; % spawn index
    T_spawn = (tspan(end) - tspan(1))/M_spawn; % time between spawn events
    tstart = tspan(1);
    tout = tstart; % output times
    yout = y0.'; % output solutions
    teout = []; % event time log
    yeout = []; % event solution log
    ieout = []; % event index log
    while tstart < tspan(end)
        % Solve until the first terminal event.
        [t,y,te,ye,ie] = ode45(@fp_derivs, tspan, y0, options);

        % Accumulate output. Dont return event times in output
        nt = length(t);
        
        if isempty(ie) || te(end) == tspan(end) % stopped due to max time      
            tout = [tout; t(2:nt)];
            yout = [yout; y(2:nt,:)];
        else % stopped due to event
            tout = [tout; t(2:nt-1)];
            yout = [yout; y(2:nt-1,:)];
        end
            
        teout = [teout; te];        % Events at tstart are never reported.
        yeout = [yeout; ye];
        ieout = [ieout; ie];

        % Check which event happened: piston or spawn?
        if ie == Np+1 % spawn time triggered event           
            
            contour = reshape(y(nt,:)', [Np, 2, Nc]);

            % insert point at piston position with correct rapidity
            contour = [xp_func(t(nt)), K0+vp_func(t(nt)); contour];

            % delete last point to conserve number of points
            contour(Np+1,:) = [];

            % set modified contour as new starting point for ODE and
            % increment spawn index
            y0 = contour(:);
            i_spawn = i_spawn + 1;
        
        else % piston touched point in upper edge
            
            y(nt,ie)        = xp_func(t(nt)); % set position equal to piston (for safety, shouldnt otherwise be necessary)
            y(nt,Np+ie)     = y(nt,Np+ie) + vp_func(t(nt)); % boost fermi edge by piston velocity
            y0              = y(nt,:);
        end

        
        % A good guess of a valid first timestep is the length of the last valid
        % timestep, so use it for faster computation.  'refine' is 4 by default.
        refine = 4;
        if length(t) <= refine
            dt_init = t(nt)-t(nt-1);
        else
            dt_init = t(nt)-t(nt-refine);
        end
        
        options = odeset(options,'InitialStep',dt_init,...
          'MaxStep',t(nt)-t(1));

        tstart = t(nt);
        tspan = [tstart tspan(tspan > tstart)];
    end


    % Reshape solution into usual contour structure
    contour_t = cell(1, length(tout));
    for i = 1:length(tout)
        ytmp        = yout(i,:)';
        ytmp        = ytmp(~isnan(ytmp)); % do not return dummy NaN
        contour     = reshape(ytmp, [], 2, Nc);        

        % ensure that points in back align with piston edge
        if contour(1,1) > xp_func(tout(i))
            contour = [ xp_func(tout(i)), K0 + vp_func(tout(i));
                        contour ];
        end

        if contour(end,2) > K0 && contour(end,1) > xp_func(tout(i)) % 2 Fermi sea state
            contour = [ contour; 
                        xp_func(tout(i)), -K0 + vp_func(tout(i)) ];
        end

        contour_t{i} = contour;

    end

        
    
    function dydt = fp_derivs(t, y)
        % Calculate velocities of contour points

        % Create auxillary countour points, representing the lower edge, 
        % which for the piston problem is constant and implicit.
        G           = reshape(y, [Np, 2, Nc]);
        G           = G(1:end-(M_spawn-i_spawn+1), :);  % remove NaN padding

        if G(end,2) > K0 % 2 Fermi sea state
            G_temp  = [ G;
                        xp_func(t), K0;
                        1e10, K0; 
                        1e10, -K0; 
                        -1e10, -K0; 
                        -1e10, K0+vp_func(t) ];
        else
            G_temp  = [ G;
                        1e10, K0; 
                        1e10, -K0; 
                        -1e10, -K0; 
                        -1e10, K0+vp_func(t) ];
        end
        
        % Calculate effective velocity
        veff        = veff_func(G_temp);
        veff        = veff(1:size(G,1),:,:); % discard auxillary points
        veff        = [ veff; 
                        nan(M_spawn-i_spawn+1, 2) ]; % add padding
            
        dydt        = veff(:);
    end


    function [value,isterminal,direction] = events(t,y)
    % Locate the time when distance between a contour point and piston
    % passes through zero in a decreasing direction and stop integration.
        dist_to_piston  = y(1:Np) - xp_func(t); % detect distance = 0
        time_to_spawn   = i_spawn*T_spawn - t; 
        
        % stop if either: time to spawn point OR piston catches up to point
        value           = [ dist_to_piston; time_to_spawn ]; 

        isterminal = ones(Np+1, 1);   % stop the integration
        direction = -ones(Np+1, 1);   % negative direction
    end
end