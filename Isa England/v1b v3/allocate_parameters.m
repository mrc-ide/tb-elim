function [p,r] = allocate_parameters(x,p,r,xi)


r.beta         = x(xi.beta);
p.betadec      = x(xi.betadec);
% p.betadec      = 0;
r.gamma_2015   = x(xi.gamma(1));
r.gamma_2020   = x(xi.gamma(2));
% r.TPT2020rec   = x(xi.r_TPT2020rec);
r.progression  = r.progression0*[1 x(xi.p_relrate) 1]*x(xi.rf_progression);
r.reactivation = r.reactivation0*[1 x(xi.p_relrate) 1]*x(xi.rf_reactivation);
r.migr         = x(xi.r_migr);
p.LTBI_in_migr = x(xi.p_LTBI_in_migr);

% nat hist rates

r.LTBI_stabil   = x(xi.LTBI_stabil);
r.Tx            = x(xi.Tx);
r.default       = x(xi.default);
r.self_cure     = x(xi.self_cure);
r.relapse       = x(xi.relapse(1):xi.relapse(3));  
r.muTB         = x(xi.muTB);
p.imm          = x(xi.imm);
