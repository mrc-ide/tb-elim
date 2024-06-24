clear all;

states   = {'U','Lf','Ls','Pf','Ps','I','I2','Tx','Rlo','Rhi','R'};
gps.born = {'dom','for'};

[i, s, d, lim] = get_addresses({states, gps.born}, [], [], [], 0);
d = char(d);

s.allI = [s.I, s.I2];

% Include the auxiliaries
names = {'inc','incsources','mort'};
lgths = [    3,           5,     1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(2,i.nstates); 
tmp(1,s.allI) = 1;
tmp(2,intersect(s.allI,s.dom)) = 1;
tmp(3,intersect(s.allI,s.for)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

% From recent infection
tmp = zeros(i.nstates);
tmp(s.allI,s.Lf) = 1;
sel.Lf2I = tmp - diag(diag(tmp));

% From recent infection, TPT
tmp = zeros(i.nstates);
tmp(s.allI,s.Pf) = 1;
sel.Pf2I = tmp - diag(diag(tmp));

% From remote infection
tmp = zeros(i.nstates);
tmp(s.allI,s.Ls) = 1;
sel.Ls2I = tmp - diag(diag(tmp));

% From remote infection, TPT
tmp = zeros(i.nstates);
tmp(s.allI,s.Ps) = 1;
sel.Ps2I = tmp - diag(diag(tmp));

% From relapse
tmp = zeros(i.nstates);
tmp(s.allI,[s.R, s.Rlo, s.Rhi]) = 1;
sel.R2I = tmp - diag(diag(tmp));

% -- Natural history parameters -------------------------------------------
r.progression  = 0.0826;
r.LTBI_stabil  = 0.872;
r.reactivation = 0.0006;

r.Tx            = 2;
r.default       = 0.01;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
% r.relapse      = [0 0 0];
% r.mu           = 1/66;                                                   % natural mortality
r.muTB         = 1/6;                                                      % TB related mortality
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
p.TPTeff       = 0.6;                                                      % Effectiveness of TPT
r.TPT          = [0 0];
r.ACF          = [0 0];
r.ACF2         = [0 0];

% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','betadec','gamma', 'p_birth','p_kLf'};                       % THESE ARE PARAMS - GAMMA IS HEALTH SYSTEM
lgths =      [1,        1,      1,        1,      1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)    = [0 40];
bds(xi.betadec,:) = [0 0.15];
bds(xi.gamma,:)   = [0 50];
bds(xi.p_birth,:) = [0 1];
bds(xi.p_kLf,:)   = [1 200];
prm.bounds = bds';


ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------

data.incd2010   = [14.1 14.6 15.1];
data.incd2020   = [6.5 7 7.5];
data.mort       = [0.28 0.3 0.32];
data.p_migrTB   = [0.708 0.728 0.748];
data.p_migrpopn = [0.138 0.168 0.198];
data.p_LTBI     = [0.15 0.2 0.25];

show = 1;
f1a = get_distribution_fns(data.incd2010, 'lognorm', show);
f1b = get_distribution_fns(data.incd2020, 'lognorm', show);
f2  = get_distribution_fns(data.mort, 'lognorm', show);
f3  = get_distribution_fns(data.p_migrTB, 'beta', show);
f4  = get_distribution_fns(data.p_migrpopn, 'beta', show);
f5  = get_distribution_fns(data.p_LTBI, 'beta', show);

% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);

save Model_setup;

