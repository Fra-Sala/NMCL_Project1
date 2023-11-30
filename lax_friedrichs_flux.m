function F = lax_friedrichs_flux(flux_phys, u, v, delta_x, k)

% Lax-Friedrichs flux: u and v have size (2,1), and correspond to qj and qj+1
F = 0.5 * (flux_phys(u) + flux_phys(v) - delta_x / k * (v - u));
    
return
