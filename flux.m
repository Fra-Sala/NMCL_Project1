function F = flux(m,h)

g = 1;

F = [m;
    m.^2/h + 1/2*g*h.^2];

return