clear
close all
clc

% implementation of the finite difference solution of the shallow-water equations

% Definition of parameters
g = 1;
% Spatial domain

xspan = [0,2];
% Temporal domain
tspan = [0,2];
h0 = @(x) 1+0.5*sin(pi*x);
m0 = @(x) 0.5*sin(pi*x);
% Number of grid points
N = 100;
% Velocity
u = 0.25;
% Number of time steps
CFL = 0.5;
k = CFL*(xspan(2)-xspan(1))/(N+1)*1/(u+sqrt(g*1.5));
K = round((tspan(end)-tspan(1))/k);

% Source function
S = @(x,t) [pi/2*(u-1)*cos(pi*(x-t)); 
            pi/2*cos(pi*(x-t))*(-u+u^2+g*h0(x-t))];

bc = 'periodic';
[h, m, tvec, xvec] = lax_friedrich(xspan, tspan, N, K, h0, m0, @flux_phys, S, bc);


for i=1:length(tvec)
    subplot(2,1,1)
    plot(xvec, h(:,i))
    title(['h(x,t) at t = ', num2str(tvec(i))])
    xlabel('x')
    ylabel('h(x,t)')
    grid on
    drawnow
    subplot(2,1,2)
    plot(xvec, m(:,i))
    title(['m(x,t) at t = ', num2str(tvec(i))])
    xlabel('x')
    ylabel('m(x,t)')
    grid on
    drawnow
end


