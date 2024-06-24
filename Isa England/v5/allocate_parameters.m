function [p,r] = allocate_parameters(x,p,r,xi)

r.beta         = x(xi.beta);
p.betadec      = x(xi.betadec);
% p.betadec      = 0;
r.gamma_2015   = x(xi.gamma(1));
r.gamma_2020   = x(xi.gamma(2));
% r.TPT2020rec   = x(xi.r_TPT2020rec);
r.progression  = r.progression0*[1 x(xi.p_relrate) 1];
r.reactivation = r.reactivation0*[1 x(xi.p_relrate) 1];
r.migr         = x(xi.r_migr);
p.LTBI_in_migr = x(xi.p_LTBI_in_migr);