clear all; load Model_setup;

xx = [5.1953    1.0520   10.1081    0.0053];
opts = odeset('NonNegative',1:i.nstates);

obj = @(x) get_objective2(x, ref, prm, gps, lhd);
[out,aux] = obj(xx);
init = aux.soln(end,:);
incdend = aux.incd(1);

% Simulate baseline
[p0,r0] = allocate_parameters(xx,p,r,xi);
M0 = make_model(p0,r0,i,s,gps);
M1 = make_model(p0,r0,i,s,gps);

geq = @(t,in) goveqs_scaleup(t, in, i, M0, M1, [2022 2025], agg, sel, r, p);
[t, soln] = ode15s(geq,[2022:2036],init);
incd(:,1) = diff(soln(:,i.aux.inc(1)),[],1)*1e5;

figure; plot(incd)
yl = ylim; yl(1) = 0; ylim(yl);