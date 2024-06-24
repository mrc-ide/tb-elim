function out = goveqs_basis2(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% New infections
lam = M.lam*invec/sum(invec);

% Full model
allmat = M.lin + lam*M.nlin;
out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

allmorts = sum(morts(:));
births = p.birth*allmorts;
out(i.U.dom) = out(i.U.dom) + births;

vec = invec(s.dom);
vec([2,3]) = vec([2,3])*p.kLf;
vec = vec/sum(vec)*(1-p.birth)*allmorts;
out(s.for) = out(s.for) + vec;

% % Migration
% out(1:i.nstates) = out(1:i.nstates) + M.migration;

% Auxiliaries
out(i.aux.inc)  = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.mort) = sum(morts(:,2));