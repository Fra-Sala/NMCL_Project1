function f = flux_phys(q)

% Physical flux function for the shallow water equations

g = 1;
h = q(1);
m = q(2);

f = [m;
    m.^2/h + 1/2*g*h.^2];

return