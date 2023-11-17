clear
close all
clc

% implementation of the finite difference solution of the shallow-water equations

% Definition of parameters
g = 1;
% Spatial domain
xa = 0;
xb = 1;
% Temporal domain
tspan = [0,1];

[h, m] = lax_friedrich(xa,xb,tspan)

