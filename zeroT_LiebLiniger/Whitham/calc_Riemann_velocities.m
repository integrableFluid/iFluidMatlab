function v = calc_Riemann_velocities(R)
    % =====================================================================
    % Purpose : Calculate characteristic velocities of Riemann invariants.
    % Input :   r       -- matrix with columns 1...4, being Riemann
    %                       invariants R1, ..., R4
    % Output:   v       -- matrix with velocities
    % =====================================================================

    M = (R(2,:)-R(1,:)).*(R(4,:)-R(3,:))./(R(4,:)-R(2,:))./(R(3,:)-R(1,:));
    
    [K, E] = ellipke(M);
    
    % general expressions
    v1 = 0.5*sum(R,1) - (R(4,:)-R(1,:)).*(R(2,:)-R(1,:)).*K./( (R(4,:)-R(1,:)).*K - (R(4,:)-R(2,:)).*E );
    v2 = 0.5*sum(R,1) + (R(3,:)-R(2,:)).*(R(2,:)-R(1,:)).*K./( (R(3,:)-R(2,:)).*K - (R(3,:)-R(1,:)).*E );
    v3 = 0.5*sum(R,1) - (R(4,:)-R(3,:)).*(R(3,:)-R(2,:)).*K./( (R(3,:)-R(2,:)).*K - (R(4,:)-R(2,:)).*E );
    v4 = 0.5*sum(R,1) + (R(4,:)-R(3,:)).*(R(4,:)-R(1,:)).*K./( (R(4,:)-R(1,:)).*K - (R(3,:)-R(1,:)).*E );
   
    v = [v1; v2; v3; v4];
    
    % limiting expressions for M = 1
    v1_M1 = 1.5*R(1,:) + 0.5*R(4,:);
    v2_M1 = 0.5*R(1,:) + 0.5*R(4,:) + R(2,:);
    v3_M1 = 0.5*R(1,:) + 0.5*R(4,:) + R(2,:);
    v4_M1 = 1.5*R(4,:) + 0.5*R(1,:);
    
    v_M1 = [v1_M1; v2_M1; v3_M1; v4_M1]; 
    
    
    % limiting expressions for M = 0
    v1_M0 = 1.5*R(1,:) + 0.5*R(2,:);
    v2_M0 = 1.5*R(2,:) + 0.5*R(1,:);
    v3_M0 = 2*R(4,:) + (R(2,:) - R(1,:)).^2 ./(2*(R(1,:) + R(2,:) - 2*R(4,:)));
    v4_M0 = 2*R(4,:) + (R(2,:) - R(1,:)).^2 ./(2*(R(1,:) + R(2,:) - 2*R(4,:)));
    
    v_M0 = [v1_M0; v2_M0; v3_M0; v4_M0]; 

    
    tol = eps;
    i_M1 = abs(M - 1) < tol; % indices for M == 1 within double precision
    i_M0 = abs(M - 0) < tol; % indices for M == 0 within double precision
    
    % replace limiting cases
    v(:,i_M1) = v_M1(:,i_M1);
    v(:,i_M0) = v_M0(:,i_M0);
    
end