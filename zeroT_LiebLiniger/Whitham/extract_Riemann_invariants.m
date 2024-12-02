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

        % find index of lower and upper DSW edge and lower contour edge
        i_low = find( contours{1,i}(:,2) < 0, 1); 
        i_up = find( diff(contours{1,i}(1:i_low-2,1))<=0 ,1);
        i_dw = i_up + find( diff(contours{1,i}(i_up:end,1))>=0 ,1) -1; 
        
        i_DSW(:,i) = [i_dw; i_up; i_low];
        
        if isempty(i_up)
            continue
        end

        z_dw = contours{1,i}(i_dw,1);
        z_up = contours{1,i}(i_up,1);

        % create spatial grid of DSW region
        z_DSW(:,i) = linspace(z_dw, z_up, M_DSW);

        % interpolate to positinal grid
        for j = 1:N1
            r1(:,i,j) = interp1(contours{j,i}(i_low:end,1), contours{j,i}(i_low:end,2), z_DSW(:,i), 'makima')';
            r2(:,i,j) = interp1(contours{j,i}(i_dw:i_low-1,1), contours{j,i}(i_dw:i_low-1,2), z_DSW(:,i), 'makima')';
            r3(:,i,j) = interp1(contours{j,i}(i_up:i_dw,1), contours{j,i}(i_up:i_dw,2), z_DSW(:,i), 'makima')';
            r4(:,i,j) = interp1(contours{j,i}(1:i_up,1), contours{j,i}(1:i_up,2), z_DSW(:,i), 'makima')';
        end
    end 

    % convert from rapidity to invariants
    r1 = hbar*r1/2/m; 
    r2 = hbar*r2/2/m;
    r3 = hbar*r3/2/m;
    r4 = hbar*r4/2/m;
    
    R = {r1, r2, r3, r4};

end