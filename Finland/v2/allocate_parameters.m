function [p,r] = allocate_parameters(x,p,r,xi)

r.beta         = x(xi.beta);
p.betadec      = x(xi.betadec);
r.gamma        = x(xi.gamma);
r.ad_mort      = x(xi.ch_mort);
r.progression  = r.progression0*[1, x(xi.relrisk)];
r.reactivation = r.reactivation0*[1, x(xi.relrisk)];