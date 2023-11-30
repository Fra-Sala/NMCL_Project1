clear
close all
clc

% SECTION 1.2
% We will use periodic boundary condition option ('peri')
bc = 'peri';
animation = 'true';

xspan = [0 2];
tspan = [0 2];

%% First set of initial conditions
h01 = @(x) 1 - 0.1 * sin(pi * x);
m01 = @(x) 0;
S1 = @(x, t) [0;
    0];

% We generate a reference solution
[h1_ex, m1_ex, tvec1_ex, xvec1_ex] = conservative_scheme(xspan, tspan, 1000, 2000, h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);


% We now proceed with a less refined solution
N = 100;

% Number of time steps
CFL = 0.5;

% k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
% K = round((tspan(end) - tspan(1)) / k);
K = 2 * N;

% Solve the problem
[h1, m1, tvec1, xvec1] = conservative_scheme(xspan, tspan, N, K, h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);

% Plot of the solution at time t = 2
figure(1)

T_f = 2;
subplot(2, 1, 1)
plot(xvec1, h1(:, end), 'LineWidth', 2)
title('$h(x, t)$ at $t = 2$', 'Interpreter', 'latex')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$h(x, T)$', 'Interpreter', 'latex')
grid on
axis equal
xlim([0 2]);
set(gca, 'Fontsize', 20)
drawnow

subplot(2, 1, 2)
plot(xvec1, m1(:, end), 'LineWidth', 2)
title('$m(x, t)$ at $t = 2$', 'Interpreter', 'latex')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$m(x, T)$', 'Interpreter', 'latex')
grid on
axis equal
xlim([0 2]);
%ylim([-0.1,0.1])
set(gca, 'Fontsize', 20)
drawnow


% We visualize the animation
if animation == "true"
    figure(2)
    for i = 1 : length(tvec1)

        subplot(2, 1, 1)
        plot(xvec1, h1(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec1, m1(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec1(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

    end
end

%% Error analysis, initial condition 1

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(1:8);

% Note that we cannot solve for small values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec1 = zeros(size(N_vec));
err_m_vec1 = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);

    % k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
    % K = round((tspan(end) - tspan(1)) / k);
    K = 2 * N;

    T_f = 2;
    [h1, m1, ~, xvec1_err] = conservative_scheme(xspan, tspan, N, K, h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);
    
    % We now want to compare h1(:, end) with h1_ex(:, end), but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xvec1
    h1_interp = interp1(xvec1_err, h1(:, end), xvec1_ex);
    m1_interp = interp1(xvec1_err, m1(:, end), xvec1_ex);
    err_h_vec1(i) = norm(h1_interp' -h1_ex(:, end)); 
    err_m_vec1(i) = norm(m1_interp' - m1_ex(:, end));

end


% Plot error
figure(3)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec1, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("err", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=2\) (case 1)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on

subplot(2,1,2)
loglog(delta_x_vec, err_m_vec1, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("err", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=2\) (case 1)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)



%% Second set of initial conditions
h02 = @(x) 1 - 0.2 * sin(2 * pi * x);
m02 = @(x) 0.5;
S2 = @(x, t) [0;
    0];

% First, a refined solution as reference "exact"
[h2_ex, m2_ex, tvec2_ex, xvec2_ex] = conservative_scheme(xspan, tspan, 1000, 2000, h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);

% Solve the problem on a less refined mesh
N = 100;
K = 200;
[h2, m2, tvec2, xvec2] = conservative_scheme(xspan, tspan, N, K, h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);



% We visualize the solution
if animation == "true"
    figure(4)
    for i = 1 : length(tvec2)

        subplot(2, 1, 1)
        plot(xvec2, h2(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec2(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec2, m2(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec2(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

    end
end


%% Error analysis, initial condition 2

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(1:8);

% Note that we cannot solve for small values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec2 = zeros(size(N_vec));
err_m_vec2 = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);

    % k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
    % K = round((tspan(end) - tspan(1)) / k);
    K = 2 * N;

    T_f = 2;
    [h2, m2, ~, xvec2_err] = conservative_scheme(xspan, tspan, N, K, h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);

    % We now want to compare h1(:, end) with h1_ex(:, end), but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xvec1
    h2_interp = interp1(xvec2_err, h2(:, end), xvec2_ex);
    m2_interp = interp1(xvec2_err, m2(:, end), xvec2_ex);

    err_h_vec2(i) = norm(h2_interp' - h2_ex(:, end)); 
    err_m_vec2(i) = norm(m2_interp' - m2_ex(:, end));

end


% Plot the error
figure(5)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec2 , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("err", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=2\) (case 2)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on

subplot(2,1,2)
loglog(delta_x_vec, err_m_vec2, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("err", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=2\) (case 2)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)
