function [p,r] = allocate_parameters(x,p,r,xi)

r.beta         = x(xi.beta);
r.gamma        = x(xi.gamma);
r.ad_mort      = x(xi.ad_mort);
r.progression  = r.progression0*[1, x(xi.relrisk)];
r.reactivation = r.reactivation0*[1, x(xi.relrisk)];