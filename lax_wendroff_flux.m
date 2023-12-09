function F = lax_wendroff_flux(flux_phys, u, v, delta_x, k)
% LAX_WENDROFF_FLUX - Computes the Lax-Wendroff flux for a given 
%                     physical flux function.
%
%   F = lax_wendroff_flux(flux_phys, u, v, delta_x, k)
%
% INPUTS:
%   flux_phys   - Function handle for the physical flux function.
%   u           - State vector at position j (qj).
%   v           - State vector at position j+1 (qj+1).
%   delta_x     - Spatial grid spacing.
%   k           - Time step size.
%
% OUTPUT:
%   F           - Lax-Wendroff flux between u and v.
%
% DESCRIPTION:
%   This function calculates the Lax-Wendroff flux between two neighboring
%   states, u and v, based on the given physical flux function (flux_phys).
%  
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% December 2023

F = 0.5 * (flux_phys(u) + flux_phys(v) - k / delta_x * ...
    flux_phys_prime((u + v) / 2) * (flux_phys(v) - flux_phys(u)));

    function f_prime = flux_phys_prime(q)
        % Compute the Jacobian of the flux
        g = 1;
        h = q(1);
        m = q(2);
        f_prime = [0, 1;
            -m.^2./(h.^2) + g*h, 2 * m./h  ];
    end
end


