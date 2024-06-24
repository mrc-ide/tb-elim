clear all; load Model_setup.mat;



xx = [6.0988    1.1860    9.8189    0.0100];
obj  = @(x) get_objective2(x, ref, prm, gps, lhd);
opts = odeset('NonNegative',1:i.nstates);

[out,aux] = obj(xx);
init = aux.soln(end,:);
incdend = aux.incd(1);

% Simulate baseline
[p0,r0] = allocate_parameters(xx,p,r,xi);
M0 = make_model(p0,r0,i,s,gps);

geq = @(t,in) goveqs_basis2(t, in, i, M0, agg, sel, r, p);
[t, soln] = ode15s(geq,[2021:2036],init);
incd(:,1) = diff(soln(:,i.aux.inc(1)),[],1);

% Simulate intervention
p1 = p0; r1 = r0;
r1.TPT = [2 2];
r1.ACF = [1 1]; 
M1 = make_model(p1,r1,i,s,gps);

geq = @(t,in) goveqs_scaleup(t, in, i, M0, M1, [2022 2025], agg, sel, r, p);
[t, soln] = ode15s(geq,[2021:2036],init);
incd(:,2) = diff(soln(:,i.aux.inc(1)),[],1);

fprintf('\n');

figure; hold on;
plot(incd*1e5);
yl = ylim; yl(1) = 0; ylim(yl);
line(xlim,0.1*[1 1],'linestyle','--');

tmp1 = diff(soln(:,i.aux.incsources),[],1)*1e5;
tmp2 = tmp1(end,:)/sum(tmp1(end,:));
mat  = reshape(tmp2,2,length(tmp2)/2);