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
    
%     geq = @(t,in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
%     [t,soln] = ode15s(geq, [2022:2031], init);
%     sdiff = diff(soln,[],1);
    %incsto(:,ii,1) = sdiff(:,i.aux.inc(1))*1e5;
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    p1 = p0; r1 = r0;
    p1.migrTPT = 1;
    M1 = make_model(p1,r1,i,s,gps);
    
    p2 = p0; r2 = r0;
    p2.migrTPT = 1;
    r2.ACF = 0.69*[1 1];
    M2 = make_model(p2,r2,i,s,gps);

    p3 = p0; r3 = r0;
    p3.migrTPT = 1;
    r3.ACF = 0.69*[1 1];
    r3.TPT = 0.69*[1 0];
    M3 = make_model(p3,r3,i,s,gps);
    
    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 0.69*[1 1];
    r4.ACF = 0.69*[1 1];
    M4 = make_model(p4,r4,i,s,gps);

    p5 = p0; r5 = r0;
    p5.migrTPT = 1;
    r5.TPT  = 0.69*[1 1];
    r5.ACF  = 0.69*[1 1];
    r5.ACF2 = [12 12];
    M5 = make_model(p5,r5,i,s,gps);
    
    %models = {M0, M2, M3, M4};
    models = {M0, M2, M3, M4, M5};    
    
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, p1, [2022 2025], agg, sel, r);
        [t,soln] = ode15s(geq, [2022:2036], init);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        mrtsto(:,ii,mi) = sdiff(:,i.aux.mort)*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;
        props(ii,:,mi) = vec/sum(vec);
    end
end
fprintf('\n');


incmat = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);
mrtmat = permute(prctile(mrtsto,[2.5,50,97.5],2),[2,1,3]);


%incidence-----------------------------------------------------------------
cols = linspecer(size(incmat,3));
figure; lw = 1.5; fs = 14;

xx = [2022:2035];
for ii = 1:size(incmat,3)
   plt = incmat(:,:,ii);
   lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
   jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
   yline(1.05,'k--');
end
yl = ylim; yl(1) = 0; ylim(yl);
xlim([2022 2035])
set(gca,'fontsize',fs);

%legend(lg, 'Baseline','ACF','ACF + domestic TPT','ACF + domestic AND migrant TPT','location','SouthWest');
legend(lg, 'Baseline','ACF','ACF + domestic TPT','ACF + domestic AND migrant TPT','+ Monthly followup post TPT or Tx', 'Elimination target','location','SouthWest');
ylabel('Rate per 100,000 population');
title('Incidence')


%mortality-----------------------------------------------------------------
cols = linspecer(size(mrtmat,3));
figure; lw = 1.5; fs = 14;

xx = [2022:2035];
for ii = 1:size(mrtmat,3)
   plt = mrtmat(:,:,ii);
   lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
   jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
end
yl = ylim; yl(1) = 0; ylim(yl);
set(gca,'fontsize',fs);

%legend(lg, 'Baseline','ACF','ACF + domestic TPT','ACF + domestic AND migrant TPT','location','SouthWest');
%legend(lg, 'Baseline','ACF','ACF + domestic TPT','ACF + domestic AND migrant TPT','+ Monthly followup post TPT or Tx','location','SouthWest');
ylabel('Rate per 100,000 population');
title('Mortality')

% Show the proportions from different sources
tmp1 = prctile(props,[2.5,50,97.5],1);
tmp2 = squeeze(tmp1(2,:,end));
mm = [sum(tmp2([1,3])), sum(tmp2([2,4])), tmp2(5)];


figure; pie(mm);
labels = {'People developing TB without history of TPT or active TB treatment', 'People developing TB after TPT', 'Relapse after active TB treatment'};
legend(labels,'Location','NorthWest','Orientation','vertical');
title('Sources of incidence in 2035 with all interventions combined')