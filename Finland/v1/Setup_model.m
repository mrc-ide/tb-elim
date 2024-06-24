clear all;

states  = {'U','Lf','Ls','I','Tx','Rlo','Rhi','R'};
gps.age = {'ad','el'};

[i, s, d, lim] = get_addresses({states, gps.age}, [], [], [], 0);
d = char(d);

% Include the auxiliaries
i.aux.inc  = i.nstates + [1:3];
i.aux.mort = i.nstates + 4;
i.nx = i.aux.mort(end);

s.infectious = [s.I];
s.prevalent  = [s.infectious, s.Tx];

% Selectors for the incidence
tmp = zeros(2,i.nstates); 
tmp(1,s.I) = 1;
tmp(2,intersect(s.I,s.ad)) = 1;
tmp(3,intersect(s.I,s.el)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,[s.Ls,s.Lf]) = 1;
tmp(s.el, s.ad) = 0;
sel.inc = tmp - diag(diag(tmp));


% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

r.Tx            = 2;
r.default       = 0.01;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
% r.relapse      = [0 0 0];
% r.mu           = 1/66;                                                        % natural mortality
r.muTB         = 1/6;                                                         % TB related mortality
p.imm          = 0.8;                                                         % Reduced susceptibility conferred by previous infection



% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','gamma','relrisk','ad_mort'};
lgths =      [1,      1,        1,        1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)    = [0 15];
bds(xi.gamma,:)   = [0 6];
bds(xi.relrisk,:) = [1 25];
bds(xi.ad_mort,:) = [0, 0.01];
prm.bounds = bds';


ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify data --------------------------------------------------------

data.incd      = [4 4.7 5.4];
data.mort      = [0.54 0.55 0.56];
data.p_eldTB   = [0.62 0.64 0.66];
data.p_eldpopn = [0.21 0.23 0.25];

show = 0;
f1 = get_distribution_fns(data.incd, 'lognorm', show);
f2 = get_distribution_fns(data.mort, 'lognorm', show);
f3 = get_distribution_fns(data.p_eldTB, 'beta', show);
f4 = get_distribution_fns(data.p_eldpopn, 'beta', show);

lhd.fn = @(incd, mort, p_eldTB, p_eldpopn) f1(incd) + f2(mort) + f3(p_eldTB) + f4(p_eldpopn);

save Model_setup;

