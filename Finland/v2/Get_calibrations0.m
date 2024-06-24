clear all; load Model_setup;

nsam = 1e3;
xsam = prm.bounds(1,:) + diff(prm.bounds,1).*lhsdesign(nsam,size(prm.bounds,2));

obj  = @(x) get_objective2(x, ref, prm, gps, lhd);
nobj = @(x) -obj(x);

mk = round(nsam/25);
out = [];
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    xx = xsam(ii,:);
    [out(ii), aux] = obj(xx);
end
fprintf('\n');
mat  = sortrows([out; 1:nsam]',-1);
ord  = mat(:,2);




for ii = 1:3
    x0 = xsam(ord(ii),:);
    for jj = 1:2
        x0 = fminsearch(nobj,x0,optimset('Display','iter'));
    end
    xs = MCMC_adaptive(obj, x0, 5e4, 1, [], [], [], 1);
    
    x2 = xs(end,:);
    cov0 = cov(xs);
    [xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, x2, 5e4, 1, [], [], cov0, 1);
    fprintf('\n');
end

save calibration_res1;


return;




return;
% x1 = fminsearch(nobj,x0,optimset('Display','iter'));

% Perform MCMC






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