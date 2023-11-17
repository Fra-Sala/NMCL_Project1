clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Implementation of the finite difference solution of the shallow-water equations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1
% Definition of parameters
g = 1;
u = 0.25;

% Spatial domain
xspan = [0, 2];

% Temporal domain
tspan = [0, 2];

% Initial conditions
h0 = @(x) 1 + 0.5 * sin(pi * x);
m0 = @(x) u * h0(x);

% Number of grid points
N = 100;

% Number of time steps
CFL = 0.5;
k = CFL*(xspan(2)-xspan(1))/(N+1)*1/(u+sqrt(g*1.5));
K = round((tspan(end)-tspan(1))/k);

% Source function
S = @(x, t) [pi/2 * (u - 1) * cos(pi * (x - t)); 
            pi/2 * cos(pi * (x - t)) * (- u + u^2 + g * h0(x - t))];

% Select a boundary condition option (either 'open' or 'peri')
bc = 'peri';
[h, m, tvec, xvec] = lax_friedrichs(xspan, tspan, N, K, h0, m0, @flux_phys, S, bc);


% We visualize the solution
for i = 1 : length(tvec)

    subplot(2, 1, 1)
    plot(xvec, h(:, i), 'LineWidth', 2)
    title(['$h(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$h(x, t)$', 'Interpreter', 'latex')
    grid on
    set(gca, 'Fontsize', 20)
    drawnow

    subplot(2, 1, 2)
    plot(xvec, m(:, i), 'LineWidth', 2)
    title(['$m(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$m(x, t)$', 'Interpreter', 'latex')
    grid on
    set(gca, 'Fontsize', 20)
    drawnow

end



%% 1.2
% We now change the initial conditions and set the boundary condition
% option to 'peri'
bc = 'peri';

h01 = @(x) 1 - 0.1 * sin(pi * x);
m01 = @(x) 0;
S1 = @(x, t) [0; 
              0];

[h1, m1, tvec1, xvec1] = lax_friedrichs(xspan, tspan, N, K, h01, m01, @flux_phys, S1, bc);

% We visualize the solution
for i = 1 : length(tvec1)

    subplot(2, 1, 1)
    plot(xvec1, h1(:, i), 'LineWidth', 2)
    title(['$h(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$h(x, t)$', 'Interpreter', 'latex')
    grid on
    set(gca, 'Fontsize', 20)
    drawnow

    subplot(2, 1, 2)
    plot(xvec1, m1(:, i), 'LineWidth', 2)
    title(['$m(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$m(x, t)$', 'Interpreter', 'latex')
    grid on
    set(gca, 'Fontsize', 20)
    drawnow
    
end


