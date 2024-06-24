clear all;

states  = {'U','Lf','Ls','Pf','Ps','I','Tx','Rlo','Rhi','R'};
gps.age = {'ch','ad'};

[i, s, d, lim] = get_addresses({states, gps.age}, [], [], [], 0);
d = char(d);

% Include the auxiliaries
i.aux.inc        = i.nstates + [1:3];
i.aux.mort       = i.nstates + 4;
i.aux.incsources = i.nstates + [5:18];
i.nx = i.aux.incsources(end);

s.infectious = [s.I];
s.prevalent  = [s.infectious, s.Tx];

% Selectors for the incidence
tmp = zeros(3,i.nstates); 
tmp(1,s.I) = 1;
tmp(2,intersect(s.I,s.ch)) = 1;
tmp(3,intersect(s.I,s.ad)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,[s.Lf,s.Ls]) = 1;
tmp(s.ad, s.ch) = 0;
sel.inc = tmp - diag(diag(tmp));

% Sources of incidence
tmp = zeros(i.nstates);
tmp(s.I,s.Lf) = 1;
sel.inc_Lf = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(s.I,s.Ls) = 1;
sel.inc_Ls = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(s.I,s.Pf) = 1;
sel.inc_Pf = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(s.I,s.Ps) = 1;
sel.inc_Ps = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(s.I, s.Rlo) = 1;
sel.inc_Rlo = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.I, s.Rhi) = 1;
sel.inc_Rhi = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.I, s.R) = 1;
sel.inc_R = sparse(tmp-diag(diag(tmp)));




% -- Natural history parameters -------------------------------------------
r.progression0  = 0.072;
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
r.TPT          = [0 0];
p.TPTeff       = 0.6;
r.ACF          = [0 0];
p.crossinf     = 0.2;

% mixmat = [1 0.2; 0.2 1];

% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','betadec','gamma','relrisk','ch_mort','r_ageing','p_relbeta_ch','p_offdiag'};
lgths =      [1,        1,      1,        1,        1,         1,             1,          1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)         = [0 30];
bds(xi.betadec,:)      = [0 0.2];
bds(xi.gamma,:)        = [0 6];
bds(xi.relrisk,:)      = [0.1 25];
bds(xi.ch_mort,:)      = [0, 0.01];
bds(xi.r_ageing,:)     = [0.02 0.3];
bds(xi.p_relbeta_ch,:) = [0.2 0.8];
bds(xi.p_offdiag,:)    = [0.01 0.4];

prm.bounds = bds';


ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify data --------------------------------------------------------

data.incd2010  = [5 10 15];
data.incd2020  = [4 7 10];
data.mort      = [0.43 0.49 0.54];
data.p_adTB   = [0.75 0.82 0.90];
data.p_adpopn = [0.75 0.82 0.90];

show = 0;
% f1a = get_distribution_fns(data.incd2010, 'lognorm', show);
% f1b = get_distribution_fns(data.incd2020, 'lognorm', show);
% % f2  = get_distribution_fns(data.mort, 'lognorm', show);
% f3  = get_distribution_fns(data.p_adTB, 'beta', show);
% f4  = get_distribution_fns(data.p_adpopn, 'beta', show);

f1a = @(x) log(unifpdf(x,data.incd2010(1),data.incd2010(2)));
f1b = @(x) log(unifpdf(x,data.incd2020(1),data.incd2020(2)));
f3  = @(x) log(unifpdf(x,data.p_adTB(1),data.p_adTB(2)));
f4  = @(x) log(unifpdf(x,data.p_adpopn(1),data.p_adpopn(2)));

lhd.fn = @(incd2010, incd2020, p_adTB, p_adpopn) f1a(incd2010) + f1b(incd2020) + f3(p_adTB) + f4(p_adpopn);

save Model_setup;