function [contour_t, tout, teout] = propagate_contour_accpiston_ode45(contour, tspan, K0, T_spawn, Xp, Vp, veff_func)
    % =====================================================================
    % Purpose : Propagate Fermi contour under the action of an accelerating
    %           piston using Matlab ODE45. Piston moves from left to right.
    % Input :   contour   -- Initial contour with each row = (z, rapid)
    %           tspan     -- Matlab ODE tspan
    %           T_spawn   -- Time between spawning new contour points
    %           Xp        -- Piston position as anonynous function of time
    %           Vp        -- Piston velocity as anonynous function of time
    %           veff_func -- Anonymous function to calculate effective
    %                           velocity for contour
    % Output:   contour_t -- Cell array with contour for each time
    %           tout      -- Corresponding evolution times
    %           teout     -- Array of times of triggered events
    % =====================================================================

    Nc      = size(contour, 3); % number of contours
    Np      = size(contour, 1); % number of discretization points
        
    y0      = contour(:);
   
    options = odeset('Events',@events);

    is   = 1; % spawn index
    tstart = tspan(1);
    tout = tstart;
    yout = y0.';
    teout = [];
    yeout = [];
    ieout = [];
    while tstart < tspan(end)
        % Solve until the first terminal event.
        [t,y,te,ye,ie] = ode45(@fp_derivs, tspan, y0, options);

        % Accumulate output. Dont return event times in output
        nt = length(t);
        
        if isempty(ie) % stopped due to max time      
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
            contour = [Xp(t(nt)), K0+Vp(t(nt)); contour];

            % delete last point to conserve number of points
            contour(Np+1,:) = [];

            % set modified contour as new starting point for ODE and
            % increment spawn index
            y0 = contour(:);
            is = is + 1;
        
        else % piston touched point in upper edge
            
            y(nt,ie)        = Xp(t(nt)); % set position equal to piston (for safety, shouldnt otherwise be necessary)
            y(nt,Np+ie)     = y(nt,Np+ie) + Vp(t(nt)); % add vel to rapid
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
        cont_upper  = reshape(yout(i,:)', [Np, 2, Nc]);
        
        % include the lower edge
        xmin        = min(Xp(tout(i)));
        xmax        = max(cont_upper(:,1));
        cont_lower  = [xmax, -K0; xmin, -K0; xmin, K0 + Vp(tout(i))];
        
        contour_t{i}  = [cont_upper; cont_lower];
        
    end
        
    
    function dydt = fp_derivs(t, y)
        % Calculate velocities of contour points

        % Create lower edge of contour, which for this problem is constant
        % and implicit
        G           = reshape(y, [Np, 2, Nc]);
        G_temp      = [G; 1e10, -K0; -1e10, -K0; -1e10, K0];
        
        % Calculate effective velocity
        veff        = veff_func(G_temp, 1:size(G,1));
        veff        = permute(veff, [1 3 2]);
        
        % Include acceleration as zeros (piston is accounted via events)
        aeff        = zeros(size(veff));
        veff(:,2,:) = aeff;
            
        dydt        = veff(:);
    end


    function [value,isterminal,direction] = events(t,y)
    % Locate the time when distance between a contour point and piston
    % passes through zero in a decreasing direction and stop integration.
        dist_to_piston  = y(1:Np) - Xp(t); % detect distance = 0
        time_to_spawn   = is*T_spawn - t; 
        
        % stop if either: time to spawn point OR piston catches up to point
        value           = [ dist_to_piston; time_to_spawn ]; 

        isterminal = ones(Np+1, 1);   % stop the integration
        direction = -ones(Np+1, 1);   % negative direction
    end
end