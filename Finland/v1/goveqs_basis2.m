function out = goveqs_basis2(t, in, i, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% New infections
lam = M.lam*invec/sum(invec);

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