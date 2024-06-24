function M = make_model(p,r,i,s,gps)

m = zeros(i.nstates);

for ib = 1:length(gps.born)
    born = gps.born{ib};
    geti = @(st) i.(st).(born);
    
    U   = geti('U');
    Lf  = geti('Lf');
    Ls  = geti('Ls');
    Pf  = geti('Pf');
    Ps  = geti('Ps');
    I   = geti('I');
    I2  = geti('I2');
    Tx  = geti('Tx');
    Rlo = geti('Rlo');
    Rhi = geti('Rhi');
    R   = geti('R');

    % Progression from 'fast' latent
    source  = Lf;
    destin  = I;
    rate    = r.progression(ib);
    m(destin, source) = m(destin, source) + rate;
    
    source  = Pf;
    destin  = I2;
    rate    = r.progression(ib)*(1-p.TPTeff);
    m(destin, source) = m(destin, source) + rate;
    
    % Stabilisation of 'fast' to 'slow' latent
    source = Lf;
    destin = Ls;
    rate   = r.LTBI_stabil;
    m(destin, source) = m(destin, source) + rate;
    
    source = Pf;
    destin = Ps;
    rate   = r.LTBI_stabil;
    m(destin, source) = m(destin, source) + rate;
    
    % Reactivation of 'slow' latent
    source  = Ls;
    destin  = I;
    rate    = r.reactivation(ib);
    m(destin, source) = m(destin, source) + rate;

    source  = Ps;
    destin  = I;
    rate    = r.reactivation(ib)*(1-p.TPTeff);
    m(destin, source) = m(destin, source) + rate;
    
    % Initiation of treatment
    source  = I;
    destins = [Tx, Rhi];
    rates   = [r.gamma, r.self_cure];
    m(destins, source) = m(destins, source) + rates';

    source  = I2;
    destins = [Tx, Rhi];
    rates   = [r.gamma, r.self_cure];
    m(destins, source) = m(destins, source) + rates';
    
    % Treatment completion or interruption
    source  = Tx;
    destins = [Rlo Rhi];
    rates   = [r.Tx, r.default];
    m(destins, source) = m(destins, source) + rates';
    
    % Relapse
    sources = [Rlo Rhi R];
    destin  = I2;
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    % Stabilisation of relapse risk
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;

    % Initiation of TPT
    source = Lf;
    destin = Pf;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;
    
    source = Ls;
    destin = Ps;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;
    
    % Case-finding
    sources = [I I2];
    destin  = Tx;
    rate    = r.ACF(ib);
    m(destin, sources) = m(destin, sources) + rate;
    
    source = I2;
    destin = Tx;
    rate    = r.ACF2(ib);
    m(destin, source) = m(destin, source) + rate;
    
end

% Transition from recent to long-term migrant status (over 5-year period)
sources = s.migr_rect;
destins = s.migr_long;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + 1/5;


% % --- Ageing process
% sources = s.ad;
% destins = s.el;
% inds = sub2ind([i.nstates, i.nstates], destins, sources);
% m(inds) = m(inds) + 1/65;

M.lin = sparse(m - diag(sum(m,1)));


% --- Nonlinear component -------------------------------------------------

m = zeros(i.nstates);
for ib = 1:length(gps.born)
    born = gps.born{ib};
    susinds = intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(born));
    m(i.Lf.(born), susinds) = 1;
end
imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
m(:,imminds) = m(:,imminds)*(1-p.imm);
M.nlin = sparse(m - diag(sum(m,1)));


% --- Force of infection --------------------------------------------------

m = zeros(1,i.nstates);
m(s.allI) = r.beta;
M.lam = sparse(m);


% --- Mortality -----------------------------------------------------------

m = zeros(i.nstates,2);
m(:,1)   = 1/83;
m(s.allI,2) = r.muTB;
M.mort   = sparse(m);


% --- Mortality -----------------------------------------------------------

% m = zeros(i.nstates,1);
% m(s.for) = r.migr;
% m(intersect(s.Lf, s.for)) = r.migr*p.kLf;
% M.migration = sparse(m);