function h = plot_rapidity_distribution(rapid, Q)
    % Plot quantity Q as function of rapidity; Q must follow fluidcell 
    % index structure. Returns plot handle.
    % If Q has quasiparticle type dependence, plot Q for each type as line
    % plot with separate colors and symbols.
    % If Q is a kernel, for each qusaiparticle type combination, plot the
    % corresponding kernel.


    N       = size(Q, 1);
    Ntypes  = size(Q, 3);

    assert( N == length(rapid), "Length of rapidity grid does not match quantity!")

    if size(Q, 2) > 1
        disp('Warning! Quantity to plot has position dependence; only first position index plotted.')
    end

    if size(Q, 4) == 1
        h = make_line_plot(rapid, Q);
    elseif size(Q, 4) == N
        make_phase_space_plot(rapid, Q);
    else
        error('Incompatible size of quantity!')
    end



    function h = make_line_plot(rapid, Q)
    
    
        colors = lines;
        symbols = {'*', '^', 'd', 'o', 's', 'x', '<', '+', 'v', '.', '>', 'hexagram'};
        skip = ceil(N/8);
    
        
        plot(rapid, squeeze(double(Q)) )
        hold on
        
        for k = 1:Ntypes % plot symbols for each type
            si  = ceil(k*skip/Ntypes);
            plot(rapid(si:skip:N), double(Q(si:skip:N, 1, k)), ...
                    'LineStyle','none', ...
                    'Marker', symbols{k}, ...
                    'Color', colors(k,:), ...
                    'MarkerSize', 4, ...
                    'MarkerFaceColor', colors(k,:) )
        end

        h = gca;
    end


    function make_phase_space_plot(rapid, Q)
    
        tiledlayout(Ntypes,Ntypes, 'TileSpacing','none','Padding','none')
        
        for i = 1:Ntypes
        for j = 1:Ntypes
            Q_plot = squeeze(double( Q(:,:,i,:,j)) );
        
            nexttile
            imagesc(rapid, rapid, Q_plot);

            h(i,j) = gca;
        end
        end
        
        colormap(bluewhitered)
    
    end


end