function out = goveqs_basis2(t, in, i, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% New infections
lam = M.lam*invec/sum(invec)*(1-p.betadec)^(max((t-2010),0));

% Full model
allmat = M.lin + lam*M.nlin;
out(1:i.nstates) = allmat*invec;

% Mortality and births
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);
out(i.U.ad) = out(i.U.ad) + sum(morts(:));

% Auxiliaries
out(i.aux.inc)  = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.mort) = sum(morts(:,2));

out(i.aux.incsources(1:2))   = agg.inc([2,3],:)*(sel.inc_Lf.*allmat)*invec;
out(i.aux.incsources(3:4))   = agg.inc([2,3],:)*(sel.inc_Ls.*allmat)*invec;
out(i.aux.incsources(5:6))   = agg.inc([2,3],:)*(sel.inc_Pf.*allmat)*invec;
out(i.aux.incsources(7:8))   = agg.inc([2,3],:)*(sel.inc_Ps.*allmat)*invec;
out(i.aux.incsources(9:10))  = agg.inc([2,3],:)*(sel.inc_Rlo.*allmat)*invec;
out(i.aux.incsources(11:12)) = agg.inc([2,3],:)*(sel.inc_Rhi.*allmat)*invec;
out(i.aux.incsources(13:14)) = agg.inc([2,3],:)*(sel.inc_R.*allmat)*invec;