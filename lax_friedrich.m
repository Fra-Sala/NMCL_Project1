function [h, m, tvec, xvec] = lax_friedrich(xspan, tspan, N, K, h0, m0, flux_phys, S, bc)

% Lax-Friedrichs method for solving the shallow water equation

tvec = linspace(tspan(1), tspan(end), K+1);
xvec = linspace(xspan(1), xspan(end), N+2);
h = zeros(N+2, K+1);
m = zeros(N+2, K+1);

delta_x = (xspan(end)-xspan(1))/(N+1);
k = (tspan(end)-tspan(1))/K;

% Set the initial conditions
h(:,1) = h0(xvec);
m(:,1) = m0(xvec);

for i=2:length(tvec)

    % Set the BC for the first position
    j =1;
    source = S(xvec(j), tvec(i));
    if bc == 'periodic'
        rhs = lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j+1, i-1); m(j+1,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(end, i-1); m(end,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k);
        h(j,i) = h(j, i-1) - k/delta_x*(rhs(1)) + k*source(1);
        m(j,i) = m(j, i-1) - k/delta_x*(rhs(2))+ k*source(2);

    elseif bc == 'open'
        
        
        h(j,i) = h(j, i-1) - k/delta_x*(lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j+1, i-1); m(j+1,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k)) + k*source(1);
        m(j,i) = m(j, i-1) - k/delta_x*(lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j+1, i-1); m(j+1,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k))+ k*source(2);
    end


    for j=2:length(xvec)-1
        source = S(xvec(j), tvec(i));
        rhs = lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j+1, i-1); m(j+1,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j-1, i-1); m(j-1,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k);
        h(j,i) = h(j, i-1) - k/delta_x*(rhs(1)) + k*source(1);
        m(j,i) = m(j, i-1) - k/delta_x*(rhs(2))+ k*source(2);
    end

    j = length(xvec);
    source = S(xvec(j), tvec(i));

    if bc == 'periodic'
        rhs = lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(1, i-1); m(1,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j-1, i-1); m(j-1,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k);
        h(j,i) = h(j, i-1) - k/delta_x*(rhs(1)) + k*source(1);
        m(j,i) = m(j, i-1) - k/delta_x*(rhs(2))+ k*source(2);

    elseif bc == 'open'
        
        source = S(xvec(j), tvec(i));
        h(j,i) = h(j, i-1) - k/delta_x*(lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j-1, i-1); m(j-1,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k)) + k*source(1);
        m(j,i) = m(j, i-1) - k/delta_x*(lax_friedrich_flux(flux_phys, [h(j, i-1); m(j,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k) - ...
            lax_friedrich_flux(flux_phys, [h(j-1, i-1); m(j-1,i-1)], [h(j, i-1); m(j,i-1)], delta_x,k))+ k*source(2);

    % Set the BC for the last position
    else
        error("You must select a valid bc option")
    end
end




    % h(2:end-1,i ) = h(2:end-1,i-1) - k/h*(flux(m(2:end-1, i-1), h(2:end-1, i-1))-flux(m(1:end-2, i-1), h(1:end-2, i-1)));
    % m(2:end-1,i ) = m(2:end-1,i-1) - k/h*(flux(m(2:end-1, i-1), h(2:end-1, i-1))-flux(m(1:end-2, i-1), h(1:end-2, i-1)));

end


