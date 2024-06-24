clear all; load calibration_res.mat;

ix0 = size(xsto,1)/2;
nx  = 200;
dx  = round(ix0/nx);
xs  = xsto(ix0:dx:end,:);

for ii = 1:size(xs,1)
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = [aux.incd2010, aux.incd2020, aux.mort, aux.p_migrTB, aux.p_migrpopn, aux.p_LTBI];
    inc(:,ii)  = aux.incd;
end
sim_pct = prctile(sims,[2.5,50,97.5],1);

% Collate data
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI];
den = alldat(:,2)';

sim_plt = sim_pct./den;
dat_plt = alldat'./den;

figure; ms = 24; hold on;

md = dat_plt(2,:); hilo = diff(dat_plt,[],1);
xx = (1:length(md)) - 0.1;
plot(xx, md, 'b.', 'markersize',ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'b', 'linestyle', 'None');

md = sim_plt(2,:); hilo = diff(sim_plt,[],1);
xx = (1:length(md)) + 0.1;
plot(xx, md, 'r.', 'markersize',ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'r', 'linestyle', 'None');

yl = ylim; yl(1) = 0; ylim(yl);

% Plot incidence timeseries
figure;
plot(inc,'Color',0.5*[1 1 1]);
yl = ylim; yl(1) = 0; ylim(yl);