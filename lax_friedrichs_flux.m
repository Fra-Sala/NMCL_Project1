function F = lax_friedrichs_flux(flux_phys, u, v, delta_x, k)
% LAX_FRIEDRICHS_FLUX - Computes the Lax-Friedrichs flux for a given 
%                       physical flux function.
%
%   F = lax_friedrichs_flux(flux_phys, u, v, delta_x, k)
%
% INPUTS:
%   flux_phys   - Function handle for the physical flux function.
%   u           - State vector at position j (qj).
%   v           - State vector at position j+1 (qj+1).
%   delta_x     - Spatial grid spacing.
%   k           - Time step size.
%
% OUTPUT:
%   F           - Lax-Friedrichs flux between u and v.
%
% DESCRIPTION:
%   This function calculates the Lax-Friedrichs flux between two neighboring
%   states, u and v, based on the given physical flux function (flux_phys).
%   The Lax-Friedrichs scheme is a numerical method often used for solving
%   hyperbolic partial differential equations.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% December 2023

% Lax-Friedrichs flux: u and v have size (2,1), and correspond to qj and qj+1
F = 0.5 * (flux_phys(u) + flux_phys(v) - delta_x / k * (v - u));
    
return
