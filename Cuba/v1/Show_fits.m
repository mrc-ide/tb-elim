clear all; load calibration_res5;

ix0 = size(xsto,1)/2;
nx  = 200;
dx  = round(ix0/nx);
xs  = xsto(ix0:dx:end,:);

for ii = 1:size(xs,1)
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = [aux.incd2010, aux.incd2020, aux.mort, aux.p_adTB, aux.p_adpopn];
    inc(:,ii)  = aux.incd;
end

alldat  = [data.incd2010; data.incd2020; data.mort; data.p_adTB; data.p_adpopn]';
sim_pct = prctile(sims,[2.5,50,97.5],1)./alldat(2,:);
alldat  = alldat./alldat(2,:);

figure; hold on;
ms = 14; lw = 1.5;

xx = (1:size(sim_pct,2))-0.1;
plt  = alldat;
md   = plt(2,:);
hilo = diff(plt,[],1);
plot(xx,md,'.','markersize',ms,'Color','r');
errorbar(xx,md,hilo(1,:),hilo(2,:),'linestyle','None','linewidth',lw,'Color','r');

xx = (1:size(sim_pct,2))+0.1;
plt  = sim_pct;
md   = plt(2,:);
hilo = diff(plt,[],1);
plot(xx,md,'.','markersize',ms,'Color','b');
errorbar(xx,md,hilo(1,:),hilo(2,:),'linestyle','None','linewidth',lw,'Color','b');
yl = ylim; yl(1) = 0; ylim(yl);

% Plot incidence timeseries
figure; hold on;
plot(inc,'Color',0.5*[1 1 1]);

% Plot data for comparison
vecs = [data.incd2010; data.incd2020]';
hilo = diff(vecs,[],1);
plot(1, vecs(2,1),'r.','markersize',24);
plot(10,vecs(2,2),'r.','markersize',24);
errorbar([1, 10], vecs(2,:), hilo(1,:), hilo(2,:),'linestyle','None');
yl = ylim; yl(1) = 0; ylim(yl);

% Plot posterior densities
figure;
fnms = fieldnames(xi);
for ii = 1:13
   subplot(2,4,ii); 
   histogram(xsto(:,ii));
   title(fnms{ii});
end