function F = lax_wendroff_flux(flux_phys, u, v, delta_x, k)

% Lax-Wendroff flux: u and v have size (2,1), and correspond to qj and qj+1

F = 0.5 * (flux_phys(u) + flux_phys(v) - k / delta_x * flux_phys_prime((u + v) / 2) * (flux_phys(v) - flux_phys(u)));

    function f_prime = flux_phys_prime(q)
    
    % Derivative of the physical flux function for the shallow water equations
    
    g = 1;
    h = q(1);
    m = q(2);
    
    f_prime = [0, 1;
                -m.^2./(h.^2) + g*h, 2 * m./h  ];
    end


end


