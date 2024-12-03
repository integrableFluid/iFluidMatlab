function [phi, drho, boost] = sample_thermal_flucts_TLL_coupled(N_samples, ...      % number of samples
                                                        N_modes, ...        % number of momentum modes
                                                        z_grid, ...         % spatial grid
                                                        v_sound, ...        % sound velocity w.r.t background
                                                        T, ...              % temperature
                                                        dens0, ...          % background density
                                                        L, ...              % system length
                                                        J, ...              % tunnel coupling strength
                                                        periodic_BC)        % flag for periodic boundary conditions
    %
    % Sample thermal phase and density fluctuations following standard
    % Luttinger Liquid quantitzation.
    % Samples around a homogeneous background, where spatial variations are
    % considered coherent fluctuations, whose amplitudes are given in the
    % array coherent_amp and angles in coherent_angle.
    %
    % NOTE: does not include zero-mode
    %
    % Also returns boost (spatial derivative) of phase phi.
    
    M           = length(z_grid);

    phi         = zeros(N_samples, M, N_modes);
    drho        = zeros(N_samples, M, N_modes);
    boost       = zeros(N_samples, M, N_modes);

    % setup the mode numbers/indices
    if periodic_BC    
        modes       = -(N_modes/2):(N_modes/2);
        modes(N_modes/2+1)= []; % remove zero-mode
    else
        modes       = 1:N_modes;
    end
    
    for j = 1:N_modes
        % calculate mode functions given boundary conditions
        if periodic_BC
            k_j         = 2*pi*modes(j)/L;
            eig_func    = exp(1i*k_j*z_grid);
            eig_deriv   = 1i*k_j*exp(1i*k_j*z_grid);
        else
            k_j         = pi*modes(j)/L; % box
            eig_func    = sqrt(2)*cos(k_j*z_grid);
            eig_deriv   = sqrt(2)*k_j*sin(k_j*z_grid);
        end
        
        % single-particle energy
        E_j         = k_j^2 + 2*J; % hbar = 2m = 1 

        % mode energy spectrum
%         eps_j       = vsound*abs(k_j); % phononic branch
        eps_j       = sqrt(E_j*(E_j + v_sound^2)); % full Bogoliubov dispersion
        
        % mode population distribution
        n_j         = 1./(exp(eps_j/T) - 1); % Bose-Einstein
%         n_j         = T/eps_k; % Reighley-Jeans
        
        f_k_plus    = 1/sqrt(L)*(eps_j/E_j).^(-1/2);
        f_k_minus   = 1/sqrt(L)*(eps_j/E_j).^(1/2);
        

        % Sample the mode creation operators
        X           = randn(N_samples, 2);        
        b_j         = sqrt(n_j/2)*(X(:,1)+1i*X(:,2)); % no coherent part
        
        % calculate phase, density and boost from the operators
        phi(:,:,j)  = 1/sqrt(4*dens0)*( (-1i)*f_k_minus.*b_j.*eig_func + 1i*conj(f_k_minus.*b_j.*eig_func) );
        drho(:,:,j) = sqrt(dens0)*( f_k_plus.*b_j.*eig_func + conj(f_k_plus.*b_j.*eig_func) );
        boost(:,:,j)= 1/sqrt(4*dens0)*( (-1i)*f_k_minus.*b_j.*eig_deriv + 1i*conj(f_k_minus.*b_j.*eig_deriv) );
    end

end