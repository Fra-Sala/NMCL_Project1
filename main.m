clear
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Implementation of the finite difference solution of the shallow-water equations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SECTION 1.1
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
k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
K = round((tspan(end) - tspan(1)) / k);

% Source function
S = @(x, t) [pi/2 * (u - 1) * cos(pi * (x - t)); 
            pi/2 * cos(pi * (x - t)) * (- u + u^2 + g * h0(x - t))];

% Here we use periodic boundary condition as the option ('peri')
bc = 'peri';

% Solve the problem
[h, m, tvec, xvec, k, delta_x] = lax_friedrichs(xspan, tspan, N, K, h0, m0, @flux_phys, S, bc);


% We visualize the solution
figure(1)
for i = 1 : length(tvec)

    subplot(2, 1, 1)
    plot(xvec, h(:, i), 'LineWidth', 2)
    hold on
    plot(xvec, h0(xvec - tvec(i)), 'Linewidth', 2)
    title(['$h(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$h(x, t)$', 'Interpreter', 'latex')
    grid on
    xlim([0 2]);
    ylim([0.4 1.6]);
    hold off
    % legend('Numerical solution', 'Exact solution', 'Interpreter', 'latex')
    set(gca, 'Fontsize', 20)
    drawnow

    subplot(2, 1, 2)
    plot(xvec, m(:, i), 'LineWidth', 2)
    hold on
    plot(xvec, u * h0(xvec - tvec(i)), 'Linewidth', 2)
    title(['$m(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$m(x, t)$', 'Interpreter', 'latex')
    grid on
    xlim([0 2]);
    ylim([0.1 0.4]);
    hold off
    % legend('Numerical solution', 'Exact solution', 'Interpreter', 'latex')
    set(gca, 'Fontsize', 20)
    drawnow

end



%% SECTION 1.2
% Again, we will use periodic boundary contion option ('peri')
bc = 'peri';

% First set of initial conditions
h01 = @(x) 1 - 0.1 * sin(pi * x);
m01 = @(x) 0;
S1 = @(x, t) [0; 
              0];


% Solve the problem
[h1, m1, tvec1, xvec1] = lax_friedrichs(xspan, tspan, N, K, h01, m01, @flux_phys, S1, bc);


% We visualize the solution
figure(2)
for i = 1 : length(tvec1)

    subplot(2, 1, 1)
    plot(xvec1, h1(:, i), 'LineWidth', 2)
    title(['$h(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$h(x, t)$', 'Interpreter', 'latex')
    grid on
    axis equal
    set(gca, 'Fontsize', 20)
    drawnow

    subplot(2, 1, 2)
    plot(xvec1, m1(:, i), 'LineWidth', 2)
    title(['$m(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$m(x, t)$', 'Interpreter', 'latex')
    grid on
    axis equal
    set(gca, 'Fontsize', 20)
    drawnow
    
end


% Second set of initial conditions
h02 = @(x) 1 - 0.2 * sin(2 * pi * x);
m02 = @(x) 0.5;
S2 = @(x, t) [0; 
              0];

% Solve the problem
[h2, m2, tvec2, xvec2] = lax_friedrichs(xspan, tspan, N, K, h02, m02, @flux_phys, S2, bc);


% We visualize the solution
figure(3)
for i = 1 : length(tvec1)

    subplot(2, 1, 1)
    plot(xvec2, h2(:, i), 'LineWidth', 2)
    title(['$h(x, t)$ at $t = $', num2str(tvec2(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$h(x, t)$', 'Interpreter', 'latex')
    grid on
    axis equal
    set(gca, 'Fontsize', 20)
    drawnow

    subplot(2, 1, 2)
    plot(xvec2, m2(:, i), 'LineWidth', 2)
    title(['$m(x, t)$ at $t = $', num2str(tvec2(i))], 'Interpreter', 'latex')
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$m(x, t)$', 'Interpreter', 'latex')
    grid on
    axis equal
    set(gca, 'Fontsize', 20)
    drawnow
    
end


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
[h, m, tvec, xvec] = lax_friedrichs(xspan, tspan, N, K, h0, m0, @flux_phys, S, bc);


% We visualize the solution
figure(4)
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
