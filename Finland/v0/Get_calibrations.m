clear all; load Model_setup;

x0 = [10 1 2 0];

obj  = @(x) get_objective2(x, ref, prm, gps, lhd);
nobj = @(x) -obj(x);

x1 = fminsearch(nobj,x0,optimset('Display','iter'));

% Perform MCMC
[xsto, outsto] = MCMC_adaptive(obj, x1, 5e4, 1, [], [], 1);



x2 = xsto2(end,:);
cov0 = cov(xsto2);
[xsto2, outsto2] = MCMC_adaptive(obj, x2, 5e4, 1, [], [], cov0, 1);
fprintf('\n');


nx  = 200;
ix0 = round(size(xsto2,1)/2);
dx  = round(size(xsto2,1)/2/nx);
xs  = xsto2(ix0:dx:end,:);

mk = round(size(xs,1)/24);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ', ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sim(ii,:) = [aux.incd, aux.mort, aux.p_eldTB, aux.p_eldpopn];
end
fprintf('\n');

save calibration_res;