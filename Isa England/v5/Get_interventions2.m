clear all; % load calibration_res.mat; 
load bestguess.mat;
load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, lhd);

% ix0 = round(size(xsto,1)/2);
% dx  = round(size(xsto,1)/2/150);
% xs  = xsto(ix0:dx:end,:);

xs = x0; 
xs(xi.gamma(2)) = 0.5;

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    r0.gamma = r0.gamma_2020;
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
    r1.TPT = 1.4*[0 1 0];
    M1 = make_model(p1,r1,i,s,gps);
    
    p2 = p1; r2 = r1;
    r2.ACF = 0.69*[1 1 1];
    M2 = make_model(p2,r2,i,s,gps);
    
    p3 = p2; r3 = r2;
    p3.migrTPT = 0.75;
    M3 = make_model(p3,r3,i,s,gps);
    
    p4 = p3; r4 = r3;
    r4.TPT = 1.4*[1 1 0];
    M4 = make_model(p4,r4,i,s,gps);
% 
%     p3 = p0; r3 = r0;
%     p3.migrTPT = 1;
%     r3.ACF = 0.69*[1 1];
%     r3.TPT = 0.69*[1 0];
%     M3 = make_model(p3,r3,i,s,gps);
%     
%     p4 = p0; r4 = r0;
%     p4.migrTPT = 1;
%     r4.TPT = 0.69*[1 1];
%     r4.ACF = 0.69*[1 1];
%     M4 = make_model(p4,r4,i,s,gps);
% 
%     p5 = p0; r5 = r0;
%     p5.migrTPT = 1;
%     r5.TPT  = 0.69*[1 1];
%     r5.ACF  = 0.69*[1 1];
%     r5.ACF2 = [12 12];
%     M5 = make_model(p5,r5,i,s,gps);
    
    %models = {M0, M2, M3, M4};
    models = {M0, M1, M2, M3, M4};    
    
    for mi = 1:length(models)
        
        if mi<4
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, p1, [2022 2025], agg, sel, r0);
        else
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, p3, [2022 2025], agg, sel, r0);
        end
        [t,soln] = ode15s(geq, [2022:2036], init);
        
        endsolsto(mi,:) = soln(end,:);
        
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


% -------------------------------------------------------------------------
% --- Plot figure of incidence and mortality impacts ----------------------

ff=figure; lw = 1.5; fs = 14;
allmat = cat(4,incmat,mrtmat);

cols = linspecer(size(allmat,3));
xx = [2022:2035];

tis = {'Incidence','Mortality'};
for is = 1:2
    subplot(1,2,is); hold on;
    
    for ii = 1:size(allmat,3)
        plt = allmat(:,:,ii,is);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
    end
    yl = ylim; yl(1) = 0; ylim(yl);
    xlim([2022 2035]);
    set(gca,'fontsize',fs);
    
    title(tis{is});
end
subplot(1,2,1);
yline(1.05,'k--');
legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
ylabel('Rate per 100,000 population');



% -------------------------------------------------------------------------
% --- Plot figure of incidence components as of 2030 ----------------------

tmp1 = reshape(props(:,:,end),3,5);                                        % Dims: 1.Dom/migr_rect/migr_long, 2.Lf,Pf,Ls,Ps,R
tmp2 = [tmp1(1,:); sum([tmp1(2,:); tmp1(3,:)],1)];                         % Dims: 1.Dom/all migr, 2.Lf,Pf,Ls,Ps,R
tmp3 = [sum(tmp2(:,[1,3]),2), sum(tmp2(:,[2,4]),2), tmp2(:,end)];          % Dims: 1.Dom/all migr, 2.All L, All P, R
tmp4 = [tmp3(1,:), tmp3(2,:)];
labels = {'UK-born without treatment history', 'UK-born after TPT', 'UK-born after TB Rx', 'Migrants without treatment history', 'Migrants after TPT', 'Migrants after TB Rx'};
figure; pie(tmp4);
legend(labels,'Location','NorthWest','Orientation','vertical');
title('Sources of incidence in 2035 with all interventions combined')


