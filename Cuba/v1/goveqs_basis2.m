function out = goveqs_basis2(t, in, i, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Find populations of children and adults, for dividing lambda by
pops = M.popnum*invec;
den = pops(1)*(M.popnum(1,:)) + pops(2)*(M.popnum(2,:));
den(den==0) = 1;
% Calculate lambdas
lam = M.lam*invec./den*(1-p.betadec)^(max((t-2010),0));

% Full model
allmat = M.lin + lam(1)*M.nlin.ch + lam(2)*M.nlin.ad;
out(1:i.nstates) = allmat*invec;

% Mortality and births
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);
out(i.U.ch) = out(i.U.ch) + sum(morts(:));

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