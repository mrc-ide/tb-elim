clear all; load calibration_res;

ix0 = round(size(xsto,1)/2);
dx = round(size(xsto,1)/2/150);
xs = xsto(ix0:dx:end,:,2);

opts = odeset('NonNegative',1:i.nstates);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ', ii/mk); end
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    init = aux.soln(end,:);
    incdend(ii) = aux.incd(1);
    
    % Simulate baseline
    [p0,r0] = allocate_parameters(xx,p,r,xi);
    M0 = make_model(p0,r0,i,s,gps);
    M1 = make_model(p0,r0,i,s,gps);
    
    geq = @(t,in) goveqs_scaleup(t, in, i, M0, M1, [2022 2025], agg, sel, r, p);
    [t, soln] = ode15s(geq,[2022:2036],init);
    incd(:,ii,1) = diff(soln(:,i.aux.inc(1)),[],1);
    
    % Simulate intervention
    p1 = p0; r1 = r0;
    
    
end
fprintf('\n');

incd_pct = prctile(incd,[2.5,50,97.5],2);