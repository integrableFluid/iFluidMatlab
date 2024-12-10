function [n_DSW, u_DSW] = calc_DSW_density(R, z_DSW, t, c)
    % =====================================================================
    % Purpose : Calculate atomic density and hydrodynamic velocity of
    %           dispersive shock wave for given Riemann invariants.
    % Input :   R       -- cell-array of Riemann invariants
    %           z_DSW   -- corresponding position grid
    %           t       -- evolution time
    %           c       -- Lieb-Liniger coupling constant
    % Output:   n       -- atomic density profile of DSW
    %           u       -- hydrodynamic velocity profile of DSW
    % =====================================================================

    % constants
    hbar = 1;
    m = 0.515;
    dens_fact = m/sqrt(c)/hbar; % r*fact --> density

    Xi = z_DSW - 0.5*(R{1}+R{2}+R{3}+R{4}).*reshape(t, 1, []); % units of length
    Xi = Xi - Xi(1,:,:); % offset to have 0 at soliton edge

    
    sn_arg = real( sqrt((R{4}-R{2}).*(R{3}-R{1})) * m/hbar ).*Xi - pi; % first convert to inverse length, times Xi is dimensionless, minus pi to have correct phase

    modul = (R{2}-R{1}).*(R{4}-R{3})./(R{4}-R{2})./(R{3}-R{1});
    modul( modul < 0 ) = 0;
    modul( modul > 1 ) = 1;
    modul( isnan(modul) ) = 1;

    % calculate density
    sn_DSW = ellipj( sn_arg, modul); 
    n_DSW = ( 0.25*(R{4}-R{3}-R{2}+R{1}).^2 + (R{4}-R{3}).*(R{2}-R{1}).*sn_DSW.^2 ) * dens_fact^2;


    % calculate velocity
    C = 1/8*(-R{1}-R{2}+R{3}+R{4}).*(-R{1}+R{2}-R{3}+R{4}).*(R{1}-R{2}-R{3}+R{4});
    u_DSW = 0.5*(R{1}+R{2}+R{3}+R{4}) + C./(n_DSW/dens_fact^2);

end