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
    init(i.U.ad) = 1 - seed;
    init(i.I.ad) = seed;
    
    geq = @(t,in)goveqs_basis2(t, in, i, M, agg, sel, r, p);
    [t0, soln0] = ode15s(geq, [0:500], init, odeset('NonNegative',1:5));
    
    dsol    = diff(soln0,[],1);
    incd    = dsol(end,i.aux.inc)*1e5;
    mort    = dsol(end,i.aux.mort)*1e5;
    p_eldTB = incd(3)/incd(1);
    
    sfin    = soln0(end,:);
    p_eldpopn = sum(sfin(s.el))/sum(sfin(1:i.nstates));

%     keyboard;
    
    out  = calfn.fn(incd(1), mort, p_eldTB, p_eldpopn);
    aux.soln = soln0;
    aux.incd = incd;
    aux.mort = mort;
    aux.p_eldTB = p_eldTB;
    aux.p_eldpopn = p_eldpopn;
end