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
    Tx  = geti('Tx');
    Rlo = geti('Rlo');
    Rhi = geti('Rhi');
    R   = geti('R');

    % Progression from 'fast' latent
    sources = [Lf Pf];
    destin  = I;
    rates   = r.progression*[1 1-p.TPTeff];
    m(destin, sources) = m(destin, sources) + rates;
    
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
    sources = [Ls Ps];
    destin  = I;
    rates   = r.reactivation*[1 1-p.TPTeff];
    m(destin, sources) = m(destin, sources) + rates;

    % Placement on TPT    
    source = Lf;
    destin = Pf;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;

    source = Ls;
    destin = Ps;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;
    
    % Initiation of treatment
    source  = I;
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
    destin  = I;
    rates   = r.relapse(ib,:);
    m(destin, sources) = m(destin, sources) + rates;
    
    % Stabilisation of relapse risk
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
    
end

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
    
    susinds = intersect([s.Pf, s.Ps],s.(born));
    m(i.Pf.(born), susinds) = 1;
end
imminds = [s.Lf, s.Ls, s.Pf, s.Ps, s.Rlo, s.Rhi, s.R];
m(:,imminds) = m(:,imminds)*(1-p.imm);
M.nlin = sparse(m - diag(sum(m,1)));


% --- Force of infection --------------------------------------------------

m = zeros(1,i.nstates);
m(s.I) = r.beta;
M.lam = sparse(m);


% --- Mortality -----------------------------------------------------------

m = zeros(i.nstates,2);
m(:,1)   = 1/83;
m(s.I,2) = r.muTB;
M.mort   = sparse(m);


% --- Mortality -----------------------------------------------------------

% m = zeros(i.nstates,1);
% m(s.for) = r.migr;
% m(intersect(s.Lf, s.for)) = r.migr*p.kLf;
% M.migration = sparse(m);