classdef effective_velocity_tabulator < handle

% Class for calculating lookup tables for effective velocity; using this,
% the effective velocity can efficiently be tabulated for states with up to
% two local Fermi sea.
    
properties (Access = public)
    
    % gridded interpolant for single Fermi sea
    veff_table_1FS

    % gridded interpolant for fluctuations around single Fermi sea
    veff_table_1FS_fluct

    % gridded interpolants for double Fermi sea
    veff_table_2FS % cell array of length 4

end % end protected properties


methods (Access = public)
    
    % Constructor
    function obj = effective_velocity_tabulator()        

    end
    
    
end % end public methods


methods (Access = public)

    function set_1FS_lookuptable( obj, veff_griddedInterpolant )
        % stores griddedInterpolant for single Fermi sea (effective
        % velocity evaluated at Fermi points)
        obj.veff_table_1FS = veff_griddedInterpolant;
    end

    function set_1FS_fluct_lookuptable( obj, veff_griddedInterpolant )
        % stores griddedInterpolant for single Fermi sea (effective
        % velocity evaluated at rapidity near Fermi points)
        obj.veff_table_1FS_fluct = veff_griddedInterpolant;
    end

    function set_2FS_lookuptable( obj, veff_griddedInterpolant_cell )
        % stores griddedInterpolant for double Fermi sea (effective
        % velocity evaluated at all 4 Fermi points)

        % input as cell array, with one cell for each FP
        obj.veff_table_2FS = veff_griddedInterpolant_cell;
    end
    
    function [veff_map, veff_table] = generate_1FS_lookuptable(obj, dq_array, c, N)
        % Calculate effective velocity (at Fermi points only) for Fermi 
        % seas of width dq, using coupling strength c and N rapidity 
        % gridpoints.
        
        M       = length(dq_array);
        dq_array= sort(dq_array);
        
        veff_table = zeros(1, M);
        
        dp = @(rapid) ones(size(rapid));
        de = @(rapid) 2*rapid;
        
        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Calculating velocity table progress:');
        cpb.start();   

        for i = 1:M
            % set symmetric Fermi points (Ground state)
            q1 = -dq_array(i)/2;
            q2 = dq_array(i)/2;
            
            % calculate effective velocity 
            [grid, weights] = create_adaptive_grids([q1; q2], N, 5);
            dp_dr = apply_dressing(dp, grid, weights, c);
            de_dr = apply_dressing(de, grid, weights, c);
            
            % store result
            veff_table(i) = de_dr(end)./dp_dr(end);

            % show progress
            cpb_text = sprintf('%d/%d velocities calculated', i, M);
            cpb.setValue(i/M);
            cpb.setText(cpb_text);
        end
        fprintf('\n')
        
        % create gridded interpolants for rapid interpolation
        veff_map = griddedInterpolant(dq_array, veff_table,'cubic');
    end


    function [veff_map, veff_table] = generate_1FS_fluct_lookuptable(obj, dq_arr, dq_max, c, N)
        % Calculate effective velocity (at rapidities up to dq_max away 
        % from Fermi points) for Fermi seas of width dq, using coupling 
        % strength c and N rapidity gridpoints.
        
        M       = length(dq_arr);
        dq_arr  = sort(dq_arr);
        
        veff_table = zeros(M, N);
        
        dp = @(rapid) ones(size(rapid));
        de = @(rapid) 2*rapid;
        
        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Calculating velocity table progress:');
        cpb.start();   

        for i = 1:M
            % set symmetric Fermi points (Ground state)
            q1 = -dq_arr(i)/2;
            q2 = dq_arr(i)/2;
            
            % setup separate grids for Fermi sea and fluctuations
            [grid_FS, weights_FS] = create_adaptive_grids([q1; q2], N, 5);
            grid_fluct  = linspace( q2-dq_max, q2+dq_max, N )';

            % solve dressing equation within Fermi sea
            dp_dr_FS    = apply_dressing(dp, grid_FS, weights_FS, c);
            de_dr_FS    = apply_dressing(de, grid_FS, weights_FS, c);

            % use self-consistency relaion of dressing to get value of
            % dressed quantity outside of Fermi sea
            dp_fluct    = dp(grid_fluct);
            de_fluct    = de(grid_fluct);
            kernel      = -1/pi * c./(c^2 + (grid_fluct - grid_FS').^2);
            dp_dr       = dp_fluct - (kernel.*weights_FS')*dp_dr_FS;
            de_dr       = de_fluct - (kernel.*weights_FS')*de_dr_FS;

            % calculate effective velocity around q2
            veff_table(i,:) = de_dr./dp_dr;

            % show progress
            cpb_text = sprintf('%d/%d velocities calculated', i, M);
            cpb.setValue(i/M);
            cpb.setText(cpb_text);
        end
        fprintf('\n')
        
        % create gridded interpolants for rapid interpolation
        delta   = linspace( -dq_max, dq_max, N ); % offset from Fermi edge
        veff_map = griddedInterpolant({dq_arr, delta}, veff_table, 'cubic');
    end
    
    
    function [veff_maps, veff_table] = generate_2FS_lookuptable(obj, dq1_arr, dq2_arr, h_arr, c, N)
        % Calculate effective velocity for 2 Fermi seas of width dq1 and
        % dq2 separated by a gap of h, using coupling strength c and N 
        % rapidity gridpoints
        
        M1    = length(dq1_arr);
        M2    = length(dq2_arr);
        M3    = length(h_arr);

        dq1_arr = sort(dq1_arr);
        dq2_arr = sort(dq2_arr);
        h_arr = sort(h_arr);
        
        veff_table = zeros(M1, M2, M3, 4); 
        
        dp = @(rapid) ones(size(rapid));
        de = @(rapid) 2*rapid;
        
        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Calculating velocity table progress:');
        cpb.start();   

        for i = 1:M1
        for j = 1:M2
        for k = 1:M3
            % set symmetric Fermi points (Ground state)
            q1 = -dq1_arr(i)/2;
            q2 = dq1_arr(i)/2;

            % second Fermi sea (only at positive rapidity for now)
            q3 = q2 + h_arr(k);
            q4 = q3 + dq2_arr(j); 
            
            % calculate effective velocity 
            [grid, weights, fp_idx] = create_adaptive_grids([q1; q2; q3; q4], N, ceil(N/5));
            dp_dr = apply_dressing(dp, grid, weights, c);
            de_dr = apply_dressing(de, grid, weights, c);
            
            veff  = de_dr(fp_idx(:))./dp_dr(fp_idx(:));

            % store result
            veff_table(i,j,k,:) = reshape(veff, [1,1,1,4]);
        end
        end
        
        % show progress
        cpb_text = sprintf('%d/%d slices of 3D table calculated', i, M1);
        cpb.setValue(i/M1);
        cpb.setText(cpb_text);
        end
        fprintf('\n')
        
        % create gridded interpolants for rapid interpolation
        veff_map1 = griddedInterpolant({dq1_arr,dq2_arr,h_arr}, veff_table(:,:,:,1),'cubic');
        veff_map2 = griddedInterpolant({dq1_arr,dq2_arr,h_arr}, veff_table(:,:,:,2),'cubic');
        veff_map3 = griddedInterpolant({dq1_arr,dq2_arr,h_arr}, veff_table(:,:,:,3),'cubic');
        veff_map4 = griddedInterpolant({dq1_arr,dq2_arr,h_arr}, veff_table(:,:,:,4),'cubic');

        veff_maps = {veff_map1, veff_map2, veff_map3, veff_map4}; 
    end


    function veff = tabulate_veff_contour(obj, Gamma, c, N, idx_eval)
        % Given a contour Gamma, tabulate the effective velocity at each
        % contour point.
        % This function detects whether 1, 2, or more local Fermi seas are
        % present, then tabulates the appropriate table. 

        if nargin < 5
            % idx_eval is the indices of Gamma at which veff is evaluated
            idx_eval = 1:size(Gamma,1);
        end
        if isempty(idx_eval)
            veff = 0;
            return
        end
        
        dp = @(rapid) ones(size(rapid));
        de = @(rapid) 2*rapid;

        % get position and rapidity of target points
        x_G   = reshape(Gamma(idx_eval,1,:),[],1);
        r_G   = reshape(Gamma(idx_eval,2,:),[],1);

        % find local Fermi points at target points
        fermi_points = find_contour_crossings(x_G, Gamma);

        veff   = zeros(length(fermi_points), 1);

        for i = 1:length(fermi_points)            
            
            switch size(fermi_points{i}, 2) % number of Fermi seas
            case 0
                % no Fermi points (should maybe throw error)
                continue

            case 1
                % find the Fermi point corresponding to contour point
                idx     = double(fermi_points{i}(2) == r_G(i)) + 1;

                % lookup velocity
                veff(i) = obj.tabulate_veff_1FS(fermi_points{i}(:), idx);

            case 2
                % find the Fermi point corresponding to contour point
                idx     = find(fermi_points{i}(:) ==  r_G(i), 1);

                % lookup velocity
                veff(i) = obj.tabulate_veff_2FS(fermi_points{i}(:), idx);

            otherwise % more than 2 Fermi seas
                % error('More than 2 local Fermi seas detected!')

                % create grids and quadratures to solve dressing equation
                [grid, weights, idx_fp] = create_adaptive_grids(fermi_points{i}(:), N, round(N/10) );
                dp_dr = apply_dressing(dp, grid, weights, c);
                de_dr = apply_dressing(de, grid, weights, c);
        
                % find entry in grid correspoding to point in contour Gamma
                idx_G = find( r_G(i) == fermi_points{i}(:), 1);
        
                veff(i)    = de_dr(idx_fp(idx_G))/dp_dr(idx_fp(idx_G));
            end  

        end    

        % reshape veff to correct format (entries matching x-vals in Gamma)
        Nc      = size(Gamma, 3); % number of contours
        Np      = length(idx_eval); % number of points in each contour
        veff    = reshape(veff, [Np, 1, Nc]);
        veff(:,2,:) = zeros(size(veff));
    end

    
    function veff = tabulate_veff_1FS(obj, fermi_points, idx)
        % Tabulate effective velocity for fermi_points assuming 1 Fermi sea
        % idx indicates which Fermi point veff should correspond to

        if nargin < 3
            idx = 1:2;
        end

        assert( length(fermi_points) == 2 )
        assert( all(ismember(idx, [1, 2])) )
        assert( ~isempty(obj.veff_table_1FS) )

        q1      = fermi_points(1);
        q2      = fermi_points(2);
                       
        % lookup effective velocity of GS and add boost
        qsign   = (-1).^idx;
        veff    = qsign.*obj.veff_table_1FS( q2-q1 ) + (q2+q1);  
    end


    function veff = tabulate_veff_2FS(obj, fermi_points, idx)
        % Tabulate effective velocity for fermi_points assuming 2 Fermi sea
        % idx indicates which Fermi point veff should correspond to

        if nargin < 3
            idx = 1:4;
        end

        assert( length(fermi_points) == 4 )
        assert( all(ismember(idx, [1, 2, 3, 4])) )
        assert( ~isempty(obj.veff_table_2FS) )

        q1      = fermi_points(1);
        q2      = fermi_points(2);
        q3      = fermi_points(3);
        q4      = fermi_points(4);

        dq1     = q2 - q1;
        dq2     = q4 - q3;
        h       = q3 - q2;
        
        % lookup effective velocity and add boost
        veff    = zeros(length(idx), 1);
        for j = 1:length(idx)
            veff(j)     = obj.veff_table_2FS{idx(j)}([dq1, dq2, h]);
        end

        veff = veff + (q2+q1);
    end


    function veff = tabulate_veff_1FS_fluct(obj, fermi_points, delta, idx)
        % Tabulate effective velocity for fluctuation at rapidity delta
        % away from the Fermi points (assumed to be 1 Fermi sea).
        % idx denotes which Fermi point to evaluate fluctuation w.r.t

        % Make sure fermi_points have correct format: each row corresponds
        % to a new pair.
        % Note, this is opposite to most of the code.
        assert( size(fermi_points, 2) == 2 )
        assert( size(fermi_points, 1) == length(delta) )
        assert( size(fermi_points, 1) == length(idx) )
        assert( all(ismember(idx, [1, 2])) )


        q1      = fermi_points(:,1);
        q2      = fermi_points(:,2);      

        % lookup effective velocity of GS and add boost
        qsign   = (-1).^idx;
        veff    = qsign.*obj.veff_table_1FS_fluct([q2-q1, qsign.*delta]) + (q2+q1);

    end

    
end % end protected methods

end % end classdef