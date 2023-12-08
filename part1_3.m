clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

%% Definition of parameters

bc = 'open';
animation = "True";

h0 = @(x) 1;
m0 = @(x) -0.5 * (x < 1);
S = @(x, t) [0;
             0];

xspan = [0, 2];
tspan = [0 0.5];

% Number of points in space and time
N = 100;
K = 200;


%% Solve the problem with Lax-Friedrichs flux
[h, m, tvec, xvec] = conservative_scheme(xspan, tspan, N, K, h0,...
    m0, @lax_friedrichs_flux, @flux_phys, S, bc);

% Reference solution with very fine mesh
[h_exf, m_exf, tvec_exf, xvec_exf] = conservative_scheme(xspan, ...
    tspan, 3000, 6000, h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);

% Animation of the solution
% Visualize the solution
if animation == "True"
    figure(1);
    for i = 1 : 20 : length(tvec)

        subplot(2, 1, 1)
        plot(xvec, h(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec, m(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec(i))],  ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

    end
end


% Compute a set of numerical solutions obtained with gradually
% decreasing mesh sizes
[h1, m1, ~, xvec1] = conservative_scheme(xspan, tspan, 100, 200,...
    h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);
[h2, m2, ~, xvec2] = conservative_scheme(xspan, tspan, 250, 500,...
    h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);
[h3, m3, ~, xvec3] = conservative_scheme(xspan, tspan, 400, 800, ...
    h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);
[h4, m4, ~, xvec4] = conservative_scheme(xspan, tspan, 800, 1600,...
    h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);
[h5, m5, ~, xvec5] = conservative_scheme(xspan, tspan, 1500, 3000, ...
    h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);

figure(2)
subplot(2, 1, 1)
plot(xvec_exf, h_exf(:, end), '-b', 'LineWidth', 3)
hold on
plot(xvec1, h1(:, end), 'LineWidth', 2)
plot(xvec2, h2(:, end), 'LineWidth', 2)
plot(xvec3, h3(:, end), 'LineWidth', 2)
plot(xvec4, h4(:, end), 'LineWidth', 2)
plot(xvec5, h5(:, end), 'LineWidth', 2)
legend('Ref: $N = 3000$', '$N = 100$', '$N = 250$', '$N = 400$',...
    '$N = 800$', '$N = 1500$', 'Interpreter', 'latex', 'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$h(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)


subplot(2, 1, 2)
plot(xvec_exf, m_exf(:, end), '-b', 'LineWidth', 3)
hold on
plot(xvec1, m1(:, end), 'LineWidth', 2)
plot(xvec2, m2(:, end), 'LineWidth', 2)
plot(xvec3, m3(:, end), 'LineWidth', 2)
plot(xvec4, m4(:, end), 'LineWidth', 2)
plot(xvec5, m5(:, end), 'LineWidth', 2)
legend('Ref: $N = 3000$', '$N = 100$', '$N = 250$', '$N = 400$',...
    '$N = 800$', '$N = 1500$', 'Interpreter', 'latex', 'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$m(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)



%% Solve the problem with Lax-Wendroff flux
[h_lw, m_lw, tvec_lw, xvec_lw] = conservative_scheme(xspan, tspan, N,...
    K, h0, m0, @lax_wendroff_flux, @flux_phys, S, bc);

% Reference solution with very fine mesh
[h_exw, m_exw, tvec_exw, xvec_exw] = conservative_scheme(xspan, tspan,...
    3000, 6000, h0, m0, @lax_wendroff_flux, @flux_phys, S, bc);

% Animation of the solution
% Visualize the solution
if animation == "True"
    figure(3);
    for i = 1 : 20 : length(tvec_lw)

        subplot(2, 1, 1)
        plot(xvec_lw, h_lw(:, i), 'LineWidth', 2)
        title(['Lax-Wendroff: $h(x, t)$ at $t = $', num2str(tvec_lw(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec_lw, m_lw(:, i), 'LineWidth', 2)
        title(['Lax-Wendroff: $m(x, t)$ at $t = $', num2str(tvec_lw(i))],...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        set(gca, 'Fontsize', 20)
        drawnow

    end
end


% We now compute a set of numerical solutions obtained with gradually
% decreasing mesh sizes
[h1, m1, ~, xvec1] = conservative_scheme(xspan, tspan, 100, 200, h0, m0,...
    @lax_wendroff_flux, @flux_phys, S, bc);
[h2, m2, ~, xvec2] = conservative_scheme(xspan, tspan, 250, 500, h0, m0,...
    @lax_wendroff_flux, @flux_phys, S, bc);
[h3, m3, ~, xvec3] = conservative_scheme(xspan, tspan, 400, 800, h0, m0,...
    @lax_wendroff_flux, @flux_phys, S, bc);
[h4, m4, ~, xvec4] = conservative_scheme(xspan, tspan, 800, 1600, h0, m0,...
    @lax_wendroff_flux, @flux_phys, S, bc);
[h5, m5, ~, xvec5] = conservative_scheme(xspan, tspan, 1500, 3000, h0, m0,...
    @lax_wendroff_flux, @flux_phys, S, bc);

figure(4)
subplot(2, 1, 1)
plot(xvec_exw, h_exw(:, end), '-b', 'LineWidth', 3)
hold on
plot(xvec1, h1(:, end), 'LineWidth', 2)
plot(xvec2, h2(:, end), 'LineWidth', 2)
plot(xvec3, h3(:, end), 'LineWidth', 2)
plot(xvec4, h4(:, end), 'LineWidth', 2)
plot(xvec5, h5(:, end), 'LineWidth', 2)
legend('Ref: $N = 3000$', '$N = 100$', '$N = 250$', '$N = 400$', ...
    '$N = 800$', '$N = 1500$', 'Interpreter', 'latex', 'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$h(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)


subplot(2, 1, 2)
plot(xvec_exw, m_exw(:, end), '-b', 'LineWidth', 3)
hold on
plot(xvec1, m1(:, end), 'LineWidth', 2)
plot(xvec2, m2(:, end), 'LineWidth', 2)
plot(xvec3, m3(:, end), 'LineWidth', 2)
plot(xvec4, m4(:, end), 'LineWidth', 2)
plot(xvec5, m5(:, end), 'LineWidth', 2)
legend('Ref: $N = 3000$', '$N = 100$', '$N = 250$', '$N = 400$',...
    '$N = 800$', '$N = 1500$', 'Interpreter', 'latex', 'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$m(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)


%% Comparison between finest solutions of Lax-Friedrichs and Lax-Wendroff
figure(5)
subplot(2, 1, 1)
plot(xvec_exw, h_exw(:, end), 'LineWidth', 2)
hold on
plot(xvec_exf, h_exf(:, end), 'LineWidth', 2)
legend('Lax-Wendroff', 'Lax-Friedrichs', 'Interpreter', 'latex', ...
    'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$h(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)


subplot(2, 1, 2)
plot(xvec_exw, m_exw(:, end), 'LineWidth', 2)
hold on
plot(xvec_exf, m_exf(:, end), 'LineWidth', 2)
legend('Lax-Wendroff', 'Lax-Friedrichs', 'Interpreter', 'latex', ...
    'Location', 'best')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$m(x, t)$', 'Interpreter', 'latex')
grid on
xlim([0 2]);
set(gca, 'Fontsize', 20)