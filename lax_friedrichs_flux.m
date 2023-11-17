function F = lax_friedrichs_flux(flux_phys, u, v, delta_x, k)

% Lax-Friedrichs flux: u and v have size (2,1), and correspond to qj and qj+1
F = (flux_phys(u) + flux_phys(v) - k/delta_x * (v - u)) / 2;
    
return
