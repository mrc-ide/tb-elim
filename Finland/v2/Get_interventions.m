clear all; load calibration_res3.mat; load Model_setup.mat;

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
    r2.ACF = 3.91202*[0 1];
    M2 = make_model(p2,r2,i,s,gps); %ACF in elderly

    p3 = p0; r3 = r0;
    p3.migrTPT = 1;
    r3.ACF = 3.91202*[0 1];
    r3.TPT = 3.91202*[0 1]; %ACF in elderly and TPT in elderly
    M3 = make_model(p3,r3,i,s,gps);
    
    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.ACF = 3.91202*[1 1]; %ACF in everyone
    r4.TPT = 3.91202*[0 1];
    M4 = make_model(p4,r4,i,s,gps);

    p5 = p0; r5 = r0;
    p5.migrTPT = 1;
    r5.TPT = 3.91202*[1 1];
    r5.ACF = 3.91202*[1 1]; %ACF and TPT in everyone
    M5 = make_model(p5,r5,i,s,gps);
    
    models = {M0, M2, M3, M4, M5};
    
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t,soln] = ode15s(geq, [2022:2041], init);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;
        props(ii,:,mi) = vec/sum(vec);
    end
end
fprintf('\n');


mat = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);

cols = linspecer(size(mat,3));
figure; lw = 3; fs = 14;

set(gca, 'FontSize', fs, 'FontWeight', 'bold');
xx = [2022:2040];
for ii = 1:size(mat,3)
   plt = mat(:,:,ii);
   lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
   jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
end
yl = ylim; yl(1) = 0; ylim(yl);
yline(0.1,'k--','LineWidth', 2);


ylabel('Incidence per 100,000 population');
xlabel('Year');

legend(lg, 'Baseline incidence','ACF in the elderly','ACF & TPT in the elderly','ACF in young adults and the elderly', 'ACF and TPT in everyone', 'location','SouthWest');

% Show the proportions from different sources
%tmp1 = prctile(props,[2.5,50,97.5],1);
%tmp2 = squeeze(tmp1(2,:,end));
% Now aggregate over age groups
%tmp3 = reshape(tmp2,2,length(tmp2)/2);
%tmp4 = sum(tmp3);
% Also add all relapses
%tmp5 = [sum(tmp4(1:2)), sum(tmp4(3:4)), sum(tmp4(5:7))];
%figure; pie(tmp5);
