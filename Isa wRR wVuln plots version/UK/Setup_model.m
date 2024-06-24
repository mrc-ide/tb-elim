clear all;

states1     = {'U'};
states2     = {'Lf','Ls','Pf','Ps','I','I2','Tx','Tx2','Rlo','Rhi','R'};
gps.born    = {'dom','migr_rect','migr_long','vuln'};
gps.strains = {'ds','rr'};

[i, s, d, lim] = get_addresses({states1, gps.born}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.born, gps.strains}, i, s, d, lim);
d = char(d);

s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.I, s.I2];
% s.migrstates = [i.U.migr_rect, i.Lf.migr_rect, i.Ls.migr_rect, i.Pf.migr_rect, i.Ps.migr_rect];
s.migrstates = intersect([s.U, s.Lf, s.Ls, s.Pf, s.Ps],s.migr_rect);
s.infectious = [s.allI, intersect(s.rr,s.Tx)];

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT'};
lgths = [    4,          24,     1,     1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(3,i.nstates); 
tmp(1,s.allI) = 1;
tmp(2,intersect(s.allI,s.migr)) = 1;
tmp(3,intersect(s.allI,s.rr))   = 1;
tmp(4,intersect(s.allI,s.vuln)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
% Remove transitions due to changing drug resistance status
tmp(s.rr, s.ds) = 0;
tmp(s.ds, s.rr) = 0;
% Remove transitions due to change in vulnerability status
tmp(s.dom, s.vuln) = 0;
tmp(s.vuln, s.dom) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

% tmp = zeros(3,i.nstates);                                                  % Rows: 1.Dom all 2.Migr all, 3.Dom RR, 4.Migr RR
% tmp(1,intersect(s.allI, s.dom)) = 1;
% tmp(2,intersect(s.allI, s.migr)) = 1;
% tmp(3,intersect(intersect(s.allI, s.dom),s.rr)) = 1;
% tmp(4,intersect(intersect(s.allI, s.migr),s.rr)) = 1;
% agg.incsources = sparse(tmp);

set1 = {s.dom, s.migr, s.vuln};
set2 = {s.ds,  s.rr};

tmp  = zeros(length(set1)*length(set2),i.nstates);
row  = 1;
for is1 = 1:length(set1)
    for is2 = 1:length(set2)
        tmp(row, intersect(s.allI, intersect(set1{is1}, set2{is2}))) = 1;
        row = row+1;
    end
end
agg.incsources = sparse(tmp);

% --- Selectors for different origins 

% Untreated TB infection
tmp = zeros(i.nstates);
tmp(s.allI, [s.Lf, s.Ls]) = 1;
sel.L2I = sparse(tmp - diag(diag(tmp)));

% Treated TB infecction
tmp = zeros(i.nstates);
tmp(s.allI, [s.Pf, s.Ps]) = 1;
sel.P2I = sparse(tmp - diag(diag(tmp)));

% Relapse, non-post-treatment
tmp = zeros(i.nstates);
tmp(s.allI, [s.R, s.Rhi]) = 1;
sel.R2I = sparse(tmp - diag(diag(tmp)));

% Relapse, post-treatment
tmp = zeros(i.nstates);
tmp(s.allI, s.Rlo) = 1;
sel.T2I = sparse(tmp - diag(diag(tmp)));



% --- People starting TPT
tmp = zeros(i.nstates);
tmp([s.Pf, s.Ps],[s.Lf, s.Ls]) = 1;
sel.nTPT = tmp - diag(diag(tmp));


% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

r.Tx            = 2;
p.RRrec         = 1;
r.RR_acqu       = 0.01;
r.Tx2           = 9/12;
r.ltfu          = 0.01;
r.ltfu2         = r.Tx2*2;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
% r.relapse      = [0 0 0];
% r.mu           = 1/66;                                                   % natural mortality
r.muTB         = 1/6;                                                      % TB related mortality
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
p.TPTeff       = [0.6 0.3];                                                      % Effectiveness of TPT
r.TPT          = [0 0 0 0];                                                % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0 0];
r.ACF2         = [0 0 0 0];


% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

% names = {'beta','betadec','gamma', 'p_birth','p_kLf','r_TPT2020rec','p_relrate','r_migr'};      
% lgths =      [1,        1,      1,         1,      1,             1,          1,       1];
% names = {'beta','betadec','gamma','r_TPT2020rec','p_relrate','r_migr','p_LTBI_in_migr'};      
% lgths =      [1,        1,      2,             1,          1,       1,               1];
names = {'beta','relbeta_RR','betadec','gamma','p_relrate','r_migr','p_LTBI_in_migr','p_RR_in_migr','r_vuln','relbeta_vuln'};      
lgths =      [1,           1,        1,      2,          2,       1,               1,             1,       1,             1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)           = [0 40];
bds(xi.relbeta_RR,:)     = [1e-4 5];
bds(xi.betadec,:)        = [0 0.15];
bds(xi.gamma,:)          = repmat([1e-4 10],2,1);
bds(xi.p_relrate,:)      = repmat([0.1 20],2,1);
bds(xi.r_migr,:)         = [0 1];
bds(xi.p_LTBI_in_migr,:) = [0 0.5];
bds(xi.p_RR_in_migr,:)   = [0 0.1];
bds(xi.r_vuln,:)         = [0 2];
bds(xi.relbeta_vuln,:)   = [0.1 20];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

prm.contmat = [1 1e-3 1e-3; 1e-3 1 1e-3; 1e-3 1e-3 1];

% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------

% data.incd2010   = [14.1 14.6 15.1];
data.incd2010   = [12 14.6 17];                                            % With broader uncertainty intervals
data.incd2020   = [6.5 7 7.5];                                             
data.incdRR2020 = [0.1 0.26 0.5];                                          % Incidence of RR-TB
data.mort       = [0.28 0.3 0.32];                                         % TB mortality, 2020
data.p_migrTB   = [0.708 0.728 0.748];                                     % Proportion contribution of migrants to overall incidence
data.p_migrpopn = [0.138 0.168 0.198];                                     % Proportion of population that is migrants
data.p_LTBI     = [0.15 0.2 0.25];                                         % Proportion of migrants with LTBI
data.p_vulnpopn = [3 5 7]/100;                                             % Proportion of UK population being vulnerable
data.p_vulnTB   = [5 10 15]/100;                                           % Proportion contribution to overall incidence

data.nTPT2019   = 1.3*[0.9 1 1.1];                                         % Number of TPT initiations in 2019, per 10^5 population

show = 0;
f1a = get_distribution_fns(data.incd2010,   'lognorm', show);
f1b = get_distribution_fns(data.incd2020,   'lognorm', show);
f1c = get_distribution_fns(data.incdRR2020, 'lognorm', show);
f2  = get_distribution_fns(data.mort,       'lognorm', show);
f3  = get_distribution_fns(data.p_migrTB,   'beta',    show);
f4  = get_distribution_fns(data.p_migrpopn, 'beta',    show);
f5  = get_distribution_fns(data.p_LTBI,     'beta',    show);
f6  = get_distribution_fns(data.p_vulnpopn, 'beta',    show);
f7  = get_distribution_fns(data.p_vulnTB,   'beta', show);

% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, nTPT2019) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(nTPT2019);
lhd.fn = @(incd2010, incd2020, incdRR2020, mort, p_migrTB, p_migrpopn, p_LTBI, p_vulnpopn, p_vulnTB) f1a(incd2010) + f1b(incd2020) + f1c(incdRR2020) + f2(mort) + ...
                                                                                                     f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + ...
                                                                                                     f6(p_vulnpopn) + f7(p_vulnTB);

save Model_setup;

