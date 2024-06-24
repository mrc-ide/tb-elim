clear all; load Model_setup;

obj  = @(x) get_objective2(x, ref, prm, gps, lhd);
nobj = @(x) -obj(x);

nsam = 1e3;
xsam = repmat(prm.bounds(1,:),nsam,1) + diff(prm.bounds).*lhsdesign(nsam,size(prm.bounds,2));
mk = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    outs(ii) = obj(xsam(ii,:));
end
% Order by fit
mat  = sortrows([outs; 1:nsam]',-1);
ord  = mat(:,2);
xord = xsam(ord,:);

x0 = fminsearch(nobj,xord(1,:),optimset('Display','iter'));

% Perform MCMC
[xsto, outsto] = MCMC_adaptive(obj, x0, 5e4, 1, [], [], [], 1);

inds = find(outsto==max(outsto));
x0 = xsto(inds(1),:);

[out, aux] = obj(x0);
sfin = aux.soln(end,:);
sum(sfin(intersect(s.for,[s.Lf,s.Ls])))/sum(sfin(s.for))


cov0 = cov(xsto);
[xsto, outsto] = MCMC_adaptive(obj, x0, 1e4, 1, [], [], cov0, 1);


save calibration_res;

return;


x2 = xsto2(end,:);
cov0 = cov(xsto2);
[xsto2, outsto2] = MCMC_adaptive(obj, x2, 5e4, 1, [], [], cov0, 1);
fprintf('\n');


nx  = 200;
ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/nx);
xs  = xsto(ix0:dx:end,:);

mk = round(size(xs,1)/24);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ', ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sim(ii,:) = [aux.incd, aux.mort, aux.p_migrTB, aux.p_migrpopn, aux.p_LTBI];
end
fprintf('\n');

save calibration_res;