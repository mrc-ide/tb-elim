function out = goveqs_basis2(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Prepare population denominators
tmp = M.denvec*invec;
den = sum(tmp.*M.denvec,1)';
den(den==0) = Inf;

% New infections
lam = M.lam*(invec./den)*(1-p.betadec)^(max((t-2010),0));
% lam = M.lam*(invec./sum(invec))*(1-p.betadec)^(max((t-2010),0));

% Full model
allmat = M.lin + lam(1)*M.nlin.dom.ds + lam(2)*M.nlin.dom.rr + ...
        lam(3)*M.nlin.migr_rect.ds + lam(4)*M.nlin.migr_rect.rr + ...
        lam(3)*M.nlin.migr_long.ds + lam(4)*M.nlin.migr_long.rr + ...
        lam(5)*M.nlin.vuln.ds      + lam(6)*M.nlin.vuln.rr;
out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

% Births into UK population
dom_morts = sum(sum(morts([s.dom,s.vuln],:)));
out(i.U.dom) = out(i.U.dom) + dom_morts;

% Migration out of UK
out(s.migr) = out(s.migr) - r.migr*invec(s.migr)/sum(invec(s.migr));

% Migration into UK
inmigr = sum(sum(morts(s.migr,:))) + r.migr;
% vec = [1-p.LTBI_in_migr, (1-p.migrTPT)*p.LTBI_in_migr*[0.02 0.98], p.migrTPT*p.LTBI_in_migr*[0.02 0.98]]';
% out(s.migrstates) = out(s.migrstates) + inmigr*vec;
out(1:i.nstates) = out(1:i.nstates) + inmigr.*M.migrentries;

% % Migration
% out(1:i.nstates) = out(1:i.nstates) + M.migration;

% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;
tmp1                  = agg.incsources*((sel.L2I.*allmat)*invec);
tmp2                  = agg.incsources*((sel.P2I.*allmat)*invec);
tmp3                  = agg.incsources*((sel.R2I.*allmat)*invec);
tmp4                  = agg.incsources*((sel.T2I.*allmat)*invec);
out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);