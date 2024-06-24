function [out, aux] = get_objective2(x, ref, prm, gps, calfn)

i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

[p,r] = allocate_parameters(x, p, r, xi);

tmp1  = [prm.bounds; x];
tmp2  = diff(tmp1([1,3,2],:),[],1);
cond1 = min(tmp2(:))<0;

if cond1
    out = -Inf;
    aux = NaN;
else
    init = zeros(1,i.nx);
    seed = 1e-5;
    init(i.U.dom)       = 1 - 0.168 - seed;
    init(i.U.migr_rect) = 0.168;
    init(i.I.dom)       = seed;
    
    p0 = p; r0 = r; p0.betadec = 0;
    r0.gamma = r0.gamma_2015;
    M0 = make_model(p0, r0, i, s, gps);
    geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, agg, sel, r0, p0);
    [t0, soln0] = ode15s(geq0, [0:5e3], init, odeset('NonNegative',1:5));

    % Scale-up of TPT amongst migrants from 2015 onwards, and increase in
    % case-finding from 2010 onwards
    p1 = p0; r1 = r0; r1.TPT = [0 r.TPT2020rec 0];
    M1 = make_model(p1, r1, i, s, gps);
    
    p2 = p; r2 = r; r2.gamma = r1.gamma_2020;
    M2 = make_model(p2, r2, i, s, gps);
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, [2015 2020; 2010 2020], i, s, p2, r2, prm, sel, agg);
    [t1, soln1] = ode15s(geq1, [2010:2020], soln0(end,:), odeset('NonNegative',1:5));
    
%     allsol = [soln0; soln1(2:end,:)];
%     allt   = [t0; t1(2:end)];
    
    dsol   = diff(soln1,[],1);
    sfin   = soln1(end,:);
    
    incd2010   = dsol(t1==2010,i.aux.inc(1))*1e5;
    incd2020   = dsol(end,i.aux.inc(1))*1e5;
    incd       = dsol(end,i.aux.inc)*1e5;
    mort       = dsol(end,i.aux.mort)*1e5;
    p_migrTB   = incd(3)/incd(1);
    
%     p_LTBI     = sum(sfin(intersect(s.migr_rect,[s.Lf, s.Ls])))/sum(sfin(s.migr_rect));
    p_LTBI     = p.LTBI_in_migr;
    p_migrpopn = sum(sfin(s.migr))/sum(sfin(1:i.nstates));

    % Number initiating TPT in 2019
    n_TPT2019  = dsol(end,i.aux.nTPT)*1e5;
    
    if incd > 0.1
        out  = calfn.fn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI);
        aux.soln       = soln1;
        aux.incd       = dsol(find(t1==2010):end,i.aux.inc(1))*1e5;
        aux.incd2010   = incd2010;
        aux.incd2020   = incd2020;
        aux.mort       = mort;
        aux.p_migrTB   = p_migrTB;
        aux.p_migrpopn = p_migrpopn;
        aux.p_LTBI     = p_LTBI;
        aux.p_migrect  = sum(sfin(s.migr_rect))/sum(sfin(1:i.nstates));
        aux.nTPT       = n_TPT2019;
    else
        out = -Inf;
        aux = NaN;
    end
    
end