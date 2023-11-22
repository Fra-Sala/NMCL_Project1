clear
close all
clc

% SECTION 1.3

%% Point a: definition of parameters

bc = 'peri';
animation = "True";
h0 = @(x) 1;
m0 = @(x) -1.5 * (x < 1);
S = @(x, t) [0;
    0];
xspan = [0, 2];
% New tspan
tspan = [0 0.5];
% Number of points in space and time
N = 100;
K = 200;


%% Point b: solution on meshes of different size


% Solve the problem
[h, m, tvec, xvec] = conservative_scheme(xspan, tspan, N, K, h0, m0,@lax_friedrichs_flux, @flux_phys, S, bc);

% Animation of the solution
% We visualize the solution
% if animation == "True"
%     figure(5)
%     for i = 1 : length(tvec)
% 
%         subplot(2, 1, 1)
%         plot(xvec, h(:, i), 'LineWidth', 2)
%         title(['$h(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
%         xlabel('$x$', 'Interpreter', 'latex')
%         ylabel('$h(x, t)$', 'Interpreter', 'latex')
%         grid on
%         axis equal
%         set(gca, 'Fontsize', 20)
%         drawnow
% 
%         subplot(2, 1, 2)
%         plot(xvec, m(:, i), 'LineWidth', 2)
%         title(['$m(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
%         xlabel('$x$', 'Interpreter', 'latex')
%         ylabel('$m(x, t)$', 'Interpreter', 'latex')
%         grid on
%         axis equal
%         set(gca, 'Fontsize', 20)
%         drawnow
% 
%     end
% end

%% Point c: Lax-Wendroff

% Solve the problem using Lax-Wendroff numerical flux
[h_lw, m_lw, tvec_lw, xvec_lw] = conservative_scheme(xspan, tspan, N, K, h0, m0,@lax_wendroff_flux, @flux_phys, S, bc);

% Animation of the solution
% We visualize the solution
if animation == "True"
    figure(5)
    for i = 1 : length(tvec)

        subplot(2, 1, 1)
        plot(xvec_lw, h_lw(:, i), 'LineWidth', 2)
        title(['Lax-Wendroff: $h(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec_lw, m_lw(:, i), 'LineWidth', 2)
        title(['Lax-Wendroff: $m(x, t)$ at $t = $', num2str(tvec(i))], 'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        axis equal
        set(gca, 'Fontsize', 20)
        drawnow

    end
end