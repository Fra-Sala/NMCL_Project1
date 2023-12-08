function [h, m, tvec, xvec, k, delta_x] = conservative_scheme(xspan, ...
    tspan, N, K, h0, m0, numerical_flux, flux_phys, S, bc)



% Implementation of the conservative scheme u_j^n+1 = u_j^n - k/h*(F_j+1/2^n-F_j-1/2^n)
% By Francesco Sala and Nicolò Viscusi

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


