function [p,r] = allocate_parameters(x,p,r,xi)

r.beta    = x(xi.beta);
p.betadec = x(xi.betadec);
r.gamma   = x(xi.gamma);
p.birth   = x(xi.p_birth);
p.kLf     = x(xi.p_kLf);