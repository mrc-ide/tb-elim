clear all; load calibration_res.mat; load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, lhd);

ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/150);
xs  = xsto(ix0:dx:end,:);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    M0 = make_model(p0,r0,i,s,gps);

    % ---------------------------------------------------------------------
    % --- Model baseline
    
    geq = @(t,in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
    [t,soln] = ode15s(geq, [2022:2031], init);
    sdiff = diff(soln,[],1);
    incsto(:,ii,1) = sdiff(:,i.aux.inc(1))*1e5;
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    p1 = p0; r1 = r0;
    p1.migrTPT = 1;
    M1 = make_model(p1,r1,i,s,gps);
    
    geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, M1, p0, p1, [2022 2025], agg, sel, r);
    [t,soln] = ode15s(geq, [2022:2031], init);
    
    sdiff = diff(soln,[],1);
    incsto(:,ii,2) = sdiff(:,i.aux.inc(1))*1e5;
end
fprintf('\n');

incpct = prctile(incsto,[2.5,50,97.5],2);

figure; hold on;
plot(incpct(:,:,1),'b');
plot(incpct(:,:,2),'r');