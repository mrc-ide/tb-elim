function M = make_model(p,r,i,s,gps)

m = zeros(i.nstates);

for ia = 1:length(gps.age)
    age = gps.age{ia};
    geti = @(st) i.(st).(age);
    
    U   = geti('U');
    Lf  = geti('Lf');
    Ls  = geti('Ls');
    I   = geti('I');
    Tx  = geti('Tx');
    Rlo = geti('Rlo');
    Rhi = geti('Rhi');
    R   = geti('R');

    % Progression from 'fast' latent
    source = Lf;
    destin = I;
    rate   = r.progression(ia);
    m(destin, source) = m(destin, source) + rate;
    
    % Stabilisation of 'fast' to 'slow' latent
    source = Lf;
    destin = Ls;
    rate   = r.LTBI_stabil;
    m(destin, source) = m(destin, source) + rate;
    
    % Reactivation of 'slow' latent
    source = Ls;
    destin = I;
    rate   = r.reactivation(ia);
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
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    % Stabilisation of relapse risk
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
   
    
end

% --- Ageing process
sources = s.ad;
destins = s.el;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + 1/65;

M.lin = sparse(m - diag(sum(m,1)));


% --- Nonlinear component -------------------------------------------------

m = zeros(i.nstates);
for ia = 1:length(gps.age)
    age = gps.age{ia};
    susinds = intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(age));
    m(i.Lf.(age), susinds) = 1;
end
imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
m(:,imminds) = m(:,imminds)*(1-p.imm);
M.nlin = sparse(m - diag(sum(m,1)));


% --- Force of infection --------------------------------------------------

m = zeros(1,i.nstates);
m(s.I) = r.beta;
M.lam = sparse(m);


% --- Mortality -----------------------------------------------------------

m = zeros(i.nstates,2);
m(s.ad,1) = r.ad_mort;
m(s.el,1) = 1/17;                % Additional life expectancy after 65 years
m(s.I,2)  = r.muTB;
M.mort = sparse(m);
