clear all; load calibration_res1; load Model_setup.mat;

ix0 = round(size(xsto,1)/2);
dx = round(size(xsto,1)/2/50);
xs = xsto(ix0:dx:end,:,2);

obj  = @(x) get_objective2(x, ref, prm, gps, lhd);
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
    
    geq = @(t,in) goveqs_basis2(t, in, i, M0, agg, sel, r, p);
    [t, soln] = ode15s(geq,[2019:2036],init);
    incd(:,ii,1) = diff(soln(:,i.aux.inc(1)),[],1);
    
    % Simulate intervention
    p1 = p0; r1 = r0;
    r1.TPT = [2 2];
    r1.ACF = [1e3 1e3];
    M1 = make_model(p1,r1,i,s,gps);

    geq = @(t,in) goveqs_scaleup(t, in, i, M0, M1, [2022 2025], agg, sel, r, p);
    [t, soln] = ode15s(geq,[2019:2036],init);
    incd(:,ii,2) = diff(soln(:,i.aux.inc(1)),[],1);
    
end
fprintf('\n');

incd_pct = prctile(incd,[2.5,50,97.5],2)*1e5;
figure; hold on;
xpts = 1:size(incd_pct,1);
cols = linspecer(size(incd_pct,3));
for ii = 1:size(incd_pct,3)
    mat = incd_pct(:,:,ii)';
    plot(xpts, mat(2,:), 'Color', cols(ii,:)); hold on;
    jbfill(xpts, mat(3,:), mat(1,:), cols(ii,:), 'None', 1, 0.2); hold on;
end
line(xlim, 0.1*[1 1], 'linestyle', '--');
yl = ylim; yl(1) = 0; ylim(yl);
