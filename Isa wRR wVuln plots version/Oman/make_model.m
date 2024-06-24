function M = make_model(p,r,i,s,gps,contmat)

m = zeros(i.nstates);

for is = 1:length(gps.strains)
    strain = gps.strains{is};
    ismdr = strcmp(strain,'rr');

    for ib = 1:length(gps.born)
        born = gps.born{ib};
        geti = @(st) i.(st).(born).(strain);

        Lf  = geti('Lf');
        Ls  = geti('Ls');
        Pf  = geti('Pf');
        Ps  = geti('Ps');
        I   = geti('I');
        I2  = geti('I2');
        Tx  = geti('Tx');
        Tx2 = geti('Tx2');
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
        rate    = r.progression(ib)*(1-p.TPTeff(is));
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
        rate    = r.reactivation(ib)*(1-p.TPTeff(is));
        m(destin, source) = m(destin, source) + rate;

        % Initiation of treatment
        pSLinit = ismdr*p.RRrec;
        source  = I;
        destins =                  [Tx,             Tx2,         Rhi];
        rates   = [r.gamma*(1-pSLinit), r.gamma*pSLinit, r.self_cure];
        m(destins, source) = m(destins, source) + rates';

        source  = I2;
        destins =                  [Tx,             Tx2,         Rhi];
        rates   = [r.gamma*(1-pSLinit), r.gamma*pSLinit, r.self_cure];
        m(destins, source) = m(destins, source) + rates';

        % Treatment completion or interruption
        source  = Tx;
        destins = [Rlo Rhi];
        rates   = [r.Tx, r.ltfu];
        m(destins, source) = m(destins, source) + rates';

        % Acquisition of drug resistance while on first-line treatment
        if ~ismdr
            source = Tx;
            destin = i.Tx.(born).rr;
            rate   = r.RR_acqu;
            m(destin, source) = m(destin, source) + rate;
        end

        % Second-line treatment
        source  = Tx2;
        destins = [Rlo Rhi];
        rates   = [r.Tx2, r.ltfu2];
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
        rate    = r.ACF(ib)*(1-pSLinit);
        m(destin, sources) = m(destin, sources) + rate;

        sources = [I I2];
        destin  = Tx2;
        rate    = r.ACF(ib)*pSLinit;
        m(destin, sources) = m(destin, sources) + rate;


        source = I2;
        destin = Tx;
        rate   = r.ACF2(ib)*(1-pSLinit);
        m(destin, source) = m(destin, source) + rate;

        source = I2;
        destin = Tx2;
        rate   = r.ACF2(ib)*pSLinit;
        m(destin, source) = m(destin, source) + rate;
    end
end

% Transition from recent to long-term migrant status (over 5-year period)
sources = s.migr_rect;
destins = s.migr_long;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + 1/5;

% Transition from dom to vulnerable population
sources = s.dom;
destins = s.vuln;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + r.vuln;

% % --- Ageing process
% sources = s.ad;
% destins = s.el;
% inds = sub2ind([i.nstates, i.nstates], destins, sources);
% m(inds) = m(inds) + 1/65;

M.lin = sparse(m - diag(sum(m,1)));


% --- Nonlinear component -------------------------------------------------

for is = 1:length(gps.strains)
    strain = gps.strains{is};
    for ib = 1:length(gps.born)
        born = gps.born{ib};

        m = zeros(i.nstates);
        susinds = intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(born));
        m(i.Lf.(born).(strain), susinds) = 1;
        
        imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
        m(:,imminds) = m(:,imminds)*(1-p.imm);
        
        M.nlin.(born).(strain) = sparse(m - diag(sum(m,1)));
    end
end

% --- Fractions for different migrant entry states ------------------------

getinds = @(st1, st2) intersect(intersect(s.migr_rect, s.(st1)), s.(st2));

m = zeros(i.nstates,1);
m(i.U.migr_rect) = 1-p.LTBI_in_migr;
m(getinds('Lf','ds')) = p.LTBI_in_migr*(1-p.migrTPT)*(1-p.RR_in_migr)*0.02;
m(getinds('Lf','rr')) = p.LTBI_in_migr*(1-p.migrTPT)*p.RR_in_migr*0.02;
m(getinds('Ls','ds')) = p.LTBI_in_migr*(1-p.migrTPT)*(1-p.RR_in_migr)*0.98;
m(getinds('Ls','rr')) = p.LTBI_in_migr*(1-p.migrTPT)*p.RR_in_migr*0.98;
m(getinds('Pf','ds')) = p.LTBI_in_migr*p.migrTPT*(1-p.RR_in_migr)*0.02;
m(getinds('Pf','rr')) = p.LTBI_in_migr*p.migrTPT*p.RR_in_migr*0.02;
m(getinds('Ps','ds')) = p.LTBI_in_migr*p.migrTPT*(1-p.RR_in_migr)*0.98;
m(getinds('Ps','rr')) = p.LTBI_in_migr*p.migrTPT*p.RR_in_migr*0.98;
M.migrentries = sparse(m);


% --- Force of infection --------------------------------------------------

getinds = @(st1, st2) intersect(intersect(s.infectious, s.(st1)), s.(st2));
contmat(end,end) = contmat(end,end)*p.relbeta_vuln;

m = zeros(6,i.nstates);                                                    % Rows: 1.Dom DS 2.Dom RR 3.Migr DS 4.Migr RR 5.Vuln DS 6.Vuln RR
m(1,getinds('dom', 'ds')) = contmat(1,1);
m(1,getinds('migr','ds')) = contmat(1,2);
m(1,getinds('vuln','ds')) = contmat(1,3);

m(2,getinds('dom', 'rr')) = contmat(1,1);
m(2,getinds('migr','rr')) = contmat(1,2);
m(2,getinds('vuln','rr')) = contmat(1,3);

m(3,getinds('dom', 'ds')) = contmat(2,1);
m(3,getinds('migr','ds')) = contmat(2,2);
m(3,getinds('vuln','ds')) = contmat(2,3);

m(4,getinds('dom', 'rr')) = contmat(2,1);
m(4,getinds('migr','rr')) = contmat(2,2);
m(4,getinds('vuln','rr')) = contmat(2,3);

m(5,getinds('dom', 'ds')) = contmat(3,1);
m(5,getinds('migr','ds')) = contmat(3,2);
m(5,getinds('vuln','ds')) = contmat(3,3);

m(6,getinds('dom', 'rr')) = contmat(3,1);
m(6,getinds('migr','rr')) = contmat(3,2);
m(6,getinds('vuln','rr')) = contmat(3,3);

% Include infectiousness
m = m*r.beta;
% Discount overall infectiousness of RR
m(:,intersect(s.rr,s.infectious)) = m(:,intersect(s.rr,s.infectious))*p.relbeta_RR;

M.lam = sparse(m);

% Additional matrix to help keep track of numbers in each group 
m = zeros(3,i.nstates);
m(1,s.dom)  = 1;
m(2,s.migr) = 1;
m(3,s.vuln) = 1;
M.denvec = sparse(m);


% --- Mortality -----------------------------------------------------------

m = zeros(i.nstates,2);
m(:,1)            = 1/83;
m(s.vuln,1)       = 1/55;
m(s.infectious,2) = r.muTB;
M.mort            = sparse(m);


% --- Mortality -----------------------------------------------------------

% m = zeros(i.nstates,1);
% m(s.for) = r.migr;
% m(intersect(s.Lf, s.for)) = r.migr*p.kLf;
% M.migration = sparse(m);