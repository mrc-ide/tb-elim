clear all;

states   = {'U','Lf','Ls','Pf','Ps','I','I2','I3','Tx','Rlo','Rhi','R'};
gps.born = {'dom','for'};

[i, s, d, lim] = get_addresses({states, gps.born}, [], [], [], 0);
d = char(d);

% Include the auxiliaries
i.aux.inc        = i.nstates + [1:6];
i.aux.mort       = i.nstates + 7;
i.aux.incsources = i.nstates + [8:21];
i.nx = i.aux.incsources(end);

s.Iall       = [s.I, s.I2, s.I3];
s.infectious = [s.I, s.I2, s.I3];
s.prevalent  = [s.infectious, s.Tx];

% Selectors for the incidence
tmp = zeros(2,i.nstates); 
tmp(1, s.Iall)                  = 1;
tmp(2,intersect(s.Iall, s.dom)) = 1;
tmp(3,intersect(s.Iall, s.for)) = 1;
tmp(4,s.I)                      = 1;
tmp(5,s.I2)                     = 1;
tmp(6,s.I3)                     = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Iall, :) = 1;
% tmp(s.el, s.ad) = 0;
sel.inc = tmp - diag(diag(tmp));

% Tracking sources of incidence
tmp = zeros(i.nstates);
tmp(s.Iall, s.Lf) = 1;
sel.inc_Lf = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.Ls) = 1;
sel.inc_Ls = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.Pf) = 1;
sel.inc_Pf = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.Ps) = 1;
sel.inc_Ps = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.Rlo) = 1;
sel.inc_Rlo = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.Rhi) = 1;
sel.inc_Rhi = sparse(tmp-diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Iall, s.R) = 1;
sel.inc_R = sparse(tmp-diag(diag(tmp)));




% -- Natural history parameters -------------------------------------------
r.progression  = 0.0826;
r.LTBI_stabil  = 0.872;
r.reactivation = 0.0006;

r.Tx            = 2;
r.default       = 0.01;

r.self_cure    = 1/6;
r.relapse      = repmat([0.032 0.14 0.0015],2,1);
% r.relapse      = [0 0 0];
% r.mu           = 1/66;                                                   % natural mortality
r.muTB         = 1/6;                                                      % TB related mortality
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection
p.migrTPT      = 0;
r.TPT          = [0 0];
p.TPTeff       = 0.6;                                                      % Effectiveness of TPT
r.followup     = 0;                                                        % Rate of followup amongst previous TB, and those on TPT
r.ACF          = 0;

% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','gamma','p_birth','p_kLf'};
lgths =      [1,      1,        1,      1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)    = [0 30];
bds(xi.gamma,:)   = [0 50];
bds(xi.p_birth,:) = [0 1];
bds(xi.p_kLf,:)   = [1 200];
prm.bounds = bds';


ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify data --------------------------------------------------------

data.incd       = [5.9 6.9 7.9];
data.mort       = [0.21 0.22 0.23];
data.p_migrTB   = [0.87 0.89 0.91];
data.p_migrpopn = [0.27 0.30 0.33];
data.p_LTBI     = [0.15 0.2 0.25];

show = 0;
f1 = get_distribution_fns(data.incd, 'lognorm', show);
f2 = get_distribution_fns(data.mort, 'lognorm', show);
f3 = get_distribution_fns(data.p_migrTB, 'beta', show);
f4 = get_distribution_fns(data.p_migrpopn, 'beta', show);
f5 = get_distribution_fns(data.p_LTBI, 'beta', show);

lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);

save Model_setup;

