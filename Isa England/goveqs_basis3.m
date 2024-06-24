function out = goveqs_basis2(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% New infections
lam = M.lam*invec/sum(invec)*(1-p.betadec)^(max((t-2010),0));

% Full model
allmat = M.lin + lam*M.nlin;
out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

% Births into UK population
dom_morts = sum(sum(morts(s.dom,:)));
out(i.U.dom) = out(i.U.dom) + dom_morts;

% Migration out of UK
out(s.migr) = out(s.migr) - r.migr*invec(s.migr)/sum(invec(s.migr));

% Migration into UK
inmigr = sum(sum(morts(s.migr,:))) + r.migr;
vec = [1-p.LTBI_in_migr, (1-p.migrTPT)*p.LTBI_in_migr*[0.02 0.98], p.migrTPT*p.LTBI_in_migr*[0.02 0.98]]';
out(s.migrstates) = out(s.migrstates) + inmigr*vec;

% % Migration
% out(1:i.nstates) = out(1:i.nstates) + M.migration;

% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;
tmp1                  = agg.incsources*((sel.Lf2I.*allmat)*invec);
tmp2                  = agg.incsources*((sel.Pf2I.*allmat)*invec);
tmp3                  = agg.incsources*((sel.Ls2I.*allmat)*invec);
tmp4                  = agg.incsources*((sel.Ps2I.*allmat)*invec);
tmp5                  = agg.incsources*((sel.R2I.*allmat)*invec);
out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4; tmp5];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);