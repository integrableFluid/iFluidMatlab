function [R, z_DSW] = extract_Riemann_invariants(contours, M_DSW)
    % =====================================================================
    % Purpose : Given Fermi contour(s), extract Riemann invariants of DSW 
    %           region with resolution specified by M_DSW. It is assumed 
    %           that the DSW is at positive rapidity.
    % Input :   contours-- cell-array of Fermi contours
    %           M_DSW   -- number of gridpoints in Riemann invariant
    % Output:   R       -- cell-array with Riemann invariants
    %           z_DSW   -- positions correpoding to each point in Riemann 
    % =====================================================================

    if ~iscell(contours)
        contours = {contours};
    end

    N1 = size(contours, 1); 
    N2 = size(contours, 2); 

    % constants
    m = 0.5;
    hbar = 1;
    
    r1 = zeros(M_DSW, N2, N1);
    r2 = zeros(M_DSW, N2, N1);
    r3 = zeros(M_DSW, N2, N1);
    r4 = zeros(M_DSW, N2, N1);
    z_DSW = zeros(M_DSW, N2);
    i_DSW = nan(3, N2);
    

    for i = 1:N2

        % Find first and last index of lower contour edge (is assumed to be
        % at negative rapidity) 
        i_low1  = find( contours{1,i}(:,2) < 0, 1, 'first'); 
        i_low2  = find( contours{1,i}(:,2) < 0, 1, 'last'); 

        % Find index of harmonic and solitonic edge of DSW 
        i_he = find( diff(contours{1,i}(1:i_low1-2,1))<=0 ,1);
        i_se = i_he + find( diff(contours{1,i}(i_he:end,1))>=0 ,1) -1; 
        
        i_DSW(:,i) = [i_se; i_he; i_low1];
        
        if isempty(i_he)
            continue
        end

        z_dw = contours{1,i}(i_se,1);
        z_up = contours{1,i}(i_he,1);

        % create spatial grid of DSW region
        z_DSW(:,i) = linspace(z_dw, z_up, M_DSW);

        % interpolate to positinal grid
        for j = 1:N1
            r1(:,i,j) = interp1(contours{j,i}(i_low1:i_low2,1), contours{j,i}(i_low1:i_low2,2), z_DSW(:,i), 'makima')';
            r2(:,i,j) = interp1(contours{j,i}(i_se:i_low1-1,1), contours{j,i}(i_se:i_low1-1,2), z_DSW(:,i), 'makima')';
            r3(:,i,j) = interp1(contours{j,i}(i_he:i_se,1), contours{j,i}(i_he:i_se,2), z_DSW(:,i), 'makima')';
            r4(:,i,j) = interp1(contours{j,i}(1:i_he,1), contours{j,i}(1:i_he,2), z_DSW(:,i), 'makima')';
        end
    end 

    % convert from rapidity to invariants
    r1 = hbar*r1/2/m; 
    r2 = hbar*r2/2/m;
    r3 = hbar*r3/2/m;
    r4 = hbar*r4/2/m;
    
    R = {r1, r2, r3, r4};

end