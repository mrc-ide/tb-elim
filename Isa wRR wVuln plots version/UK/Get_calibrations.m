clear all; load Model_setup; % load calibration_res_prev cov0;

obj  = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);
nobj = @(x) -obj(x);

nsam = 5e3; 
xsam = repmat(prm.bounds(1,:),nsam,1) + diff(prm.bounds).*lhsdesign(nsam,size(prm.bounds,2));

% obj(xsam(1,:));

% xx = xsam(1,:);
% [p,r] = allocate_parameters(xx, p, r, xi);
% init = zeros(1,i.nx);
% seed = 1e-5;
% init(i.U.dom)       = 1 - 0.168 - seed;
% init(i.U.migr_rect) = 0.168;
% init(i.I.dom)       = seed;
% 
% p0 = p; r0 = r; p0.betadec = 0;
% M0 = make_model(p0, r0, i, s, gps);
% geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, agg, sel, r0, p0);
% [t0, soln0] = ode15s(geq0, [0:1e4], init, odeset('NonNegative',1:5));

mk = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    outs(ii) = obj(xsam(ii,:));
end

% Order by fit
mat  = sortrows([outs; 1:nsam]',-1);
ord  = mat(:,2);
xord = xsam(ord,:);

load calibration_res x0_init;
x0 = x0_init;

% x0 = fminsearch(nobj,xord(1,:),optimset('Display','iter'));
x0 = fminsearch(nobj,x0,optimset('Display','iter'));
x0 = fminsearch(nobj,x0,optimset('Display','iter'));

x0_init = x0;
save calibration_isa_uk;

return;

% % Perform MCMC
% [xsto, outsto] = MCMC_adaptive(obj, x0, 5e4, 1, [], [], cov0, 1);
% 
% inds = find(outsto==max(outsto));
% x0 = xsto(inds(1),:);
% 
% [out, aux] = obj(x0);
% sfin = aux.soln(end,:);
% sum(sfin(intersect(s.migr,[s.Lf,s.Ls])))/sum(sfin(s.migr))
% 
% 
% cov0 = cov(xsto);
% [xsto, outsto] = MCMC_adaptive(obj, x0, 1e4, 1, [], [], cov0, 1);

cov0 = [];

nreps = 4;
niter = [1, 1, 1, 5]*1e4;
for ii = 1:nreps
    [xsto, outsto] = MCMC_adaptive2(obj, x0, niter(ii), 1, cov0, 1);
    inds = find(outsto==max(outsto));
    x0 = xsto(inds(1),:);
    cov0 = cov(xsto);
end




[xsto, outsto] = MCMC_adaptive(obj, x0, niter(ii), 1, [], [], cov0, 1);
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