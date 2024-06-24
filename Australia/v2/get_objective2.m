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
    M = make_model(p,r,i,s,gps);
    
    init = zeros(1,i.nx);
    seed = 1e-5;
    init(i.U.dom) = 1 - seed;
    init(i.I.dom) = seed;
    
    geq = @(t,in)goveqs_basis2(t, in, i, s, M, agg, sel, r, p);
    [t0, soln0] = ode15s(geq, [0:2020], init, odeset('NonNegative',1:5));
    
    dsol     = diff(soln0,[],1);
    incd2010 = dsol(t0==2010,i.aux.inc(1))*1e5;
    incd2020 = dsol(end,i.aux.inc(1))*1e5;
    incd     = dsol(end,i.aux.inc)*1e5;
    mort     = dsol(end,i.aux.mort)*1e5;
    p_migrTB = incd(3)/incd(1);
    
    sfin = soln0(end,:);
    p_LTBI   = sum(sfin(intersect(s.for,[s.Lf, s.Ls])))/sum(sfin(s.for));
    
    sfin    = soln0(end,:);
    p_migrpopn = sum(sfin(s.for))/sum(sfin(1:i.nstates));
    
    if incd > 0.1
        out  = calfn.fn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI);
        aux.soln       = soln0;
        aux.incd       = dsol(find(t0==2010):end,i.aux.inc(1))*1e5;
        aux.incd2010   = incd2010;
        aux.incd2020   = incd2020;
        aux.mort       = mort;
        aux.p_migrTB   = p_migrTB;
        aux.p_migrpopn = p_migrpopn;
        aux.p_LTBI     = p_LTBI;
    else
        out = -Inf;
        aux = NaN;
    end
    
end