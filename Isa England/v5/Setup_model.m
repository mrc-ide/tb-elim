clear all;

states   = {'U','Lf','Ls','Pf','Ps','I','I2','Tx','Rlo','Rhi','R'};
gps.born = {'dom','migr_rect','migr_long'};

[i, s, d, lim] = get_addresses({states, gps.born}, [], [], [], 0);
d = char(d);

s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.I, s.I2];
s.migrstates = [i.U.migr_rect, i.Lf.migr_rect, i.Ls.migr_rect, i.Pf.migr_rect, i.Ps.migr_rect];

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT'};
lgths = [    3,          15,     1,     1];
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
tmp(3,intersect(s.allI,s.migr)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

tmp = zeros(3,i.nstates);                                                  % Rows: 1.UK dom 2. UK migr, recent 3. UK migr, long-term 
tmp(1,intersect([s.I, s.I2], s.dom)) = 1;
tmp(2,intersect([s.I, s.I2], s.migr_rect)) = 1;
tmp(3,intersect([s.I, s.I2], s.migr_long)) = 1;
agg.incsources = sparse(tmp);

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


% --- People starting TPT
tmp = zeros(i.nstates);
tmp([s.Pf, s.Ps],[s.Lf, s.Ls]) = 1;
sel.nTPT = tmp - diag(diag(tmp));


% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

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
r.TPT          = [0 0 0];                                                  % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0];
r.ACF2         = [0 0 0];


% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

% names = {'beta','betadec','gamma', 'p_birth','p_kLf','r_TPT2020rec','p_relrate','r_migr'};      
% lgths =      [1,        1,      1,         1,      1,             1,          1,       1];
% names = {'beta','betadec','gamma','r_TPT2020rec','p_relrate','r_migr','p_LTBI_in_migr'};      
% lgths =      [1,        1,      2,             1,          1,       1,               1];
names = {'beta','betadec','gamma','p_relrate','r_migr','p_LTBI_in_migr'};      
lgths =      [1,        1,      2,          1,       1,               1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)           = [0 40];
bds(xi.betadec,:)        = [0 0.15];
bds(xi.gamma,:)          = repmat([0 10],2,1);
% bds(xi.p_birth,:)      = [0 1];
% bds(xi.p_kLf,:)        = [1 30];
% bds(xi.r_TPT2020rec,:)   = [0 0.01];
bds(xi.p_relrate,:)      = [1 20];
bds(xi.r_migr,:)         = [0 1];
bds(xi.p_LTBI_in_migr,:) = [0 0.5];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------

% data.incd2010   = [14.1 14.6 15.1];
data.incd2010   = [12 14.6 17];                                            % With broader uncertainty intervals
data.incd2020   = [6.5 7 7.5];                                             
data.mort       = [0.28 0.3 0.32];                                         % UK TB mortality, 2020
data.p_migrTB   = [0.708 0.728 0.748];                                     % Proportion contribution of migrants to UK incidence
data.p_migrpopn = [0.138 0.168 0.198];                                     % Proportion of UK population that is migrants
data.p_LTBI     = [0.15 0.2 0.25];                                         % Proportion of migrants with LTBI
data.nTPT2019   = 1.3*[0.9 1 1.1];                                         % Number of TPT initiations in 2019, per 10^5 population

show = 0;
f1a = get_distribution_fns(data.incd2010,   'lognorm', show);
f1b = get_distribution_fns(data.incd2020,   'lognorm', show);
f2  = get_distribution_fns(data.mort,       'lognorm', show);
f3  = get_distribution_fns(data.p_migrTB,   'beta',    show);
f4  = get_distribution_fns(data.p_migrpopn, 'beta',    show);
f5  = get_distribution_fns(data.p_LTBI,     'beta',    show);
f6  = get_distribution_fns(data.nTPT2019,   'lognorm', show);

% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, nTPT2019) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(nTPT2019);
lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);

save Model_setup;

