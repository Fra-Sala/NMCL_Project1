function [h, m, tvec, xvec, k, delta_x] = conservative_scheme(xspan, ...
    tspan, N, K, h0, m0, numerical_flux, flux_phys, S, bc)


% CONSERVATIVE_SCHEME - Implements the conservative scheme
%                       for solving hyperbolic PDEs.
%
%   [h, m, tvec, xvec, k, delta_x] = conservative_scheme(xspan, ...
%           tspan, N, K, h0, m0, numerical_flux, flux_phys, S, bc)
%
% INPUTS:
%   xspan           - Spatial domain [x_start, x_end].
%   tspan           - Temporal domain [t_start, t_end].
%   N               - Number of spatial grid points.
%   K               - Number of temporal grid points.
%   h0              - Function handle for initial water depth.
%   m0              - Function handle for initial discharge.
%   numerical_flux - Function handle for numerical flux computation.
%   flux_phys       - Function handle for the physical flux function.
%   S               - Source term function.
%   bc              - Boundary condition: 'peri' (periodic), 'open' (open).
%
% OUTPUTS:
%   h               - Water depth (matrix) over space and time.
%   m               - Discharge (matrix) over space and time.
%   tvec            - Temporal grid vector.
%   xvec            - Spatial grid vector.
%   k               - Time step size.
%   delta_x         - Spatial grid spacing.
%
% DESCRIPTION:
%   This function implements the conservative scheme for solving hyperbolic
%   partial differential equations (PDEs) that model shallow water flow. 
%   It evolves the water depth (h) and discharge (m) over a specified
%   spatial and temporal domain using the conservative scheme formula:
%
%   u_j^n+1 = u_j^n - k/h * (F_j+1/2^n - F_j-1/2^n) + k * S_j^n,
%
%   where u represents [h; m], F is the flux, S is the source term, k is
%   the time step, and delta_x is the spatial grid spacing. The boundary 
%   conditions (bc) can be set to 'peri' (periodic), 'open' (open).
%
% Authors: Francesco Sala and Nicolo' Viscusi
% December 2023


tvec    = linspace(tspan(1), tspan(end), K + 1);
xvec    = linspace(xspan(1), xspan(end), N + 1);
h       = zeros(N + 1, K + 1);
m       = zeros(N + 1, K + 1);

delta_x = (xspan(2) - xspan(1)) / N;
k       = (tspan(2) - tspan(1)) / K;

% Set the initial conditions
h(:, 1) = h0(xvec);
m(:, 1) = m0(xvec);


for i = 2 : length(tvec)

    j = 1;
    source = S(xvec(j), tvec(i));

    if bc == 'peri'

        rhs = numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(j + 1, i - 1); m(j + 1, i - 1)], delta_x, k) - ...
            numerical_flux(flux_phys, [h(end - 1, i - 1); ...
            m(end - 1, i - 1)], [h(j, i - 1); m(j, i - 1)], delta_x, k);
        h(j, i) = h(j, i - 1) - k/delta_x * rhs(1) + k * source(1);
        m(j, i) = m(j, i - 1) - k/delta_x * rhs(2) + k * source(2);

    elseif bc == 'open'

        rhs = numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(j + 1, i - 1); m(j + 1,i - 1)], delta_x, k) - ...
            numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(j, i - 1); m(j, i - 1)], delta_x, k);
        h(j, i) = h(j, i - 1) - k/delta_x * rhs(1) + k * source(1);
        m(j, i) = m(j, i - 1) - k/delta_x * rhs(2) + k * source(2);

    end

    for j = 2 : (length(xvec) - 1)

        source = S(xvec(j), tvec(i));
        rhs = numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(j + 1, i - 1); m(j + 1, i - 1)], delta_x, k) - ...
            numerical_flux(flux_phys, [h(j - 1, i - 1);  ...
            m(j - 1, i - 1)], [h(j, i - 1); m(j, i - 1)], delta_x, k);
        h(j, i) = h(j, i - 1) - k/delta_x * rhs(1) + k * source(1);
        m(j, i) = m(j, i - 1) - k/delta_x * rhs(2) + k * source(2);

    end

    j = length(xvec);
    source = S(xvec(j), tvec(i));

    if bc == 'peri'

        rhs = numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(2, i - 1); m(2, i - 1)], delta_x, k) - ...
            numerical_flux(flux_phys, [h(j - 1, i - 1); ...
            m(j - 1, i - 1)], [h(j, i - 1); m(j, i - 1)], delta_x, k);
        h(j, i) = h(j, i - 1) - k/delta_x * rhs(1) + k * source(1);
        m(j, i) = m(j, i - 1) - k/delta_x * rhs(2) + k * source(2);

    elseif bc == 'open'
            
        rhs = numerical_flux(flux_phys, [h(j, i - 1); m(j, i - 1)], ...
            [h(j, i - 1); m(j, i - 1)], delta_x, k) - ...
            numerical_flux(flux_phys, [h(j - 1, i - 1);  ...
            m(j - 1, i - 1)], [h(j, i - 1); m(j, i - 1)], delta_x, k);
        h(j, i) = h(j, i - 1) - k/delta_x * rhs(1) + k * source(1);
        m(j, i) = m(j, i - 1) - k/delta_x * rhs(2) + k * source(2);

    end

end

end


