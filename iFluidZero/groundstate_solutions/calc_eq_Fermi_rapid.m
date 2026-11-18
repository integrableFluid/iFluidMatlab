function K = calc_eq_Fermi_rapid(mu,c,N)
    % =====================================================================
    % Purpose : Calculate Fermi rapidity of ground-state Bose gas given
    %           distribution of chemical potential using Yang-Yang TBA. 
    % Input :   mu          -- Array of chemical potential values
    %           c           -- Lieb-Liniger coupling constant
    %           N           -- Number of rapidity gridpoints
    % Output:   K           -- Array of Fermi rapidities. If no Fermi
    %                          exists at a point, K(i)=NaN there.
    % =====================================================================

    K = zeros(size(mu));
    % For each mu_i, find corresponding Fermi rapidity K_i by fitting K
    % such that chemical potential computed from c,K makes mu_i.
    for i = 1:length(mu)
        % Fermi rapidity only exists if mu_i >= 0.
        if mu(i) >= 0
            K(i) = fzero(@(k) calcMu(c,k,N) - mu(i), [sqrt(mu(i)), ...
                pi*mu(i)/(8*c) + sqrt(mu(i) + (pi*mu(i)/(8*c))^2)]);
        else 
            K(i) = nan;
        end
    end
    
    % Function computing chemical potential mu0 from given values of c, K,
    % and number of rapidity points in Fermi sea
    function mu = calcMu(c, K, N_rap)
        [r, rw] = create_adaptive_grids([-K;K], N_rap, 10);
        Q_f = @(rap) [ones(size(rap)), rap.^2];
        dress = apply_dressing(Q_f, r, rw, c);
        mu = dress(1,2)/dress(1,1);
        end
end