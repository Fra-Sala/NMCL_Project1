clear
close all
clc

%% SECTION 1.3
% We will now consider discontinuous initial conditions with the open
% option ('open')
bc = 'open';

h0 = @(x) 1;
m0 = @(x) -1.5 * (x < 1);
S = @(x, t) [0;
    0];

% New tspan
tspan = [0 0.5];

% Solve the problem
[h, m, tvec, xvec] = conservative_scheme(xspan, tspan, N, K, h0, m0,@lax_friedrichs_flux, @flux_phys, S, bc);


% We visualize the solution
if animation == "true"
    figure(5)
    for i = 1 : length(tvec)

        subplot(2, 1, 1)
        plot(xvec, h(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec, m(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        set(gca, 'Fontsize', 20)
        drawnow

    end
end