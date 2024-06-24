function M = make_model(p,r,i,s,gps)

m = zeros(i.nstates);

for ia = 1:length(gps.age)
    age = gps.age{ia};
    geti = @(st) i.(st).(age);
    
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
    rates   = r.progression(ia)*[1 1-p.TPTeff];
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
    rates   = r.reactivation(ia)*[1 1-p.TPTeff];
    m(destin, sources) = m(destin, sources) + rates;

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
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    % Stabilisation of relapse risk
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
    
    % Preventive therapy
    source = Lf; 
    destin = Pf;
    rate   = r.TPT(ia);
    m(destin, source) = m(destin, source) + rate;

    source = Ls;
    destin = Ps;
    rate   = r.TPT(ia);
    m(destin, source) = m(destin, source) + rate;

    % Case-finding
    source = I;
    destin = Tx;
    rate   = r.ACF(ia);
    m(destin, source) = m(destin, source) + rate;
end

% --- Ageing process
sources = s.ch;
destins = s.ad;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + r.ageing;

M.lin = sparse(m - diag(sum(m,1)));


% --- Nonlinear component -------------------------------------------------

for ia = 1:length(gps.age)
    age = gps.age{ia};
    m = zeros(i.nstates);
    susinds = intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(age));
    m(i.Lf.(age), susinds) = 1;
    imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
    m(:,imminds) = m(:,imminds)*(1-p.imm);
    M.nlin.(age) = sparse(m - diag(sum(m,1)));
end


% --- Force of infection --------------------------------------------------

m = zeros(2,i.nstates);
m(:,s.infectious) = r.beta;
m(1,intersect(s.infectious, s.ad)) = m(1,intersect(s.infectious, s.ad))*p.offdiag;
m(2,intersect(s.infectious, s.ch)) = m(2,intersect(s.infectious, s.ch))*p.offdiag;
m(:,intersect(s.infectious,s.ch))  = m(:,intersect(s.infectious,s.ch))*p.relbeta_ch;
M.lam = sparse(m);

m = zeros(2,i.nstates);
m(1,s.ch) = 1;
m(2,s.ad) = 1;
M.popnum = sparse(m);

% --- Mortality -----------------------------------------------------------

m = zeros(i.nstates,2);
m(s.ch,1) = r.ch_mort;
m(s.ad,1) = 1/58;                % CUBA: Additional life expectancy after 15 years
m(s.I,2)  = r.muTB;
M.mort = sparse(m);
