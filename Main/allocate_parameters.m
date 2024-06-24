function [p,r] = allocate_parameters(x,p,r,xi)

r.beta         = x(xi.beta);
p.relbeta_RR   = x(xi.relbeta_RR);
p.betadec      = x(xi.betadec);
r.gamma_2015   = x(xi.gamma(1));
r.gamma_2020   = x(xi.gamma(2));
r.progression  = r.progression0*[1 x(xi.p_relrate(1)) 1 x(xi.p_relrate(2))];
r.reactivation = r.reactivation0*[1 x(xi.p_relrate(1)) 1 x(xi.p_relrate(2))];
r.migr         = x(xi.r_migr);
p.LTBI_in_migr = x(xi.p_LTBI_in_migr);
p.RR_in_migr   = x(xi.p_RR_in_migr);
r.vuln         = x(xi.r_vuln);
p.relbeta_vuln = 1;
r.ageing       = x(xi.r_ageing);
r.ch_mort      = x(xi.ch_mort);