clear all; load calib_isa.mat; load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

set1 = {'ds','rr'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
[inci, incs, incd, lim] = get_addresses({set3, set2, set1}, [], [], [], 0);


midpt = true; 
if midpt
    % inds = find(outsto==max(outsto));
    % xs = xsto(inds(1),:);
    xs = x0_init;
else
    ix0 = size(xsto,1)/2;
    nx  = 200;
    dx  = round(ix0/nx);
    xs  = xsto(ix0:dx:end,:);
end

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    r0.gamma = r0.gamma_2020;
    M0 = make_model(p0,r0,i,s,gps,prm.contmat);
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    TPTcov = -log(0.5); %TPTcov = 100;
    ACFcov = -log(0.5); %ACFcov = 100;

    % TPT in existing migrants
    ra = r0; pa = p0;
    ra.TPT = TPTcov*[0 1 1 0];
    Ma = make_model(pa,ra,i,s,gps,prm.contmat);
    
    % TPT at point of entry
    rb = ra; pb = pa;
    pb.migrTPT = 1;
    Mb = make_model(pb,rb,i,s,gps,prm.contmat);

    % TPT in contacts
    rb2 = rb; pb2 = pb;
    rb2.progression  = rb2.progression*0.9;
    rb2.reactivation = rb2.reactivation*0.9;
    Mb2 = make_model(pb2,rb2,i,s,gps,prm.contmat);

    % TPT in vulnerables
    rc = rb2; pc = pb2;
    rc.TPT = TPTcov*[0 1 1 1];
    Mc = make_model(pc,rc,i,s,gps,prm.contmat);

    % ACF in migrants and vulnerables
    rd = rc; pd = pc;
    rd.ACF  = ACFcov*[0 1 1 1];
    rd.ACF2 = rd.ACF;
    Md = make_model(pd,rd,i,s,gps,prm.contmat);

    % ACF in general population
    re = rd; pe = pd;
    re.ACF  = ACFcov*[1 1 1 1];
    re.ACF2 = rd.ACF;
    re.TPT  = TPTcov*[1 1 1 1];
    Me = make_model(pe,re,i,s,gps,prm.contmat);
    

    models = {M0, Ma, Mb, Mb2, Mc, Md, Me};    
    for mi = 1:length(models)
        
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pa, [2024 2029], agg, sel, r0);
        [t,soln] = ode15s(geq, [2022:2041], init);
        
        endsolsto(mi,:) = soln(end,:);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;

        incstoRR(:,ii,mi) = sdiff(:,i.aux.inc(3))*1e5;
        
        mat = sdiff(:,i.aux.incsources)*1e5;
        incsto2(:,ii,mi) = sum(mat(:,[incs.L,incs.R]),2);
        % incsto2(:,ii,mi) = sum(mat(:,:),2);

        mrtsto(:,ii,mi) = sdiff(:,i.aux.mort)*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;
        props(ii,:,mi) = vec/sum(vec);
        vec1 = vec;
    end
end
fprintf('\n');


incmat   = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);
incmat2  = permute(prctile(incsto2,[2.5,50,97.5],2),[2,1,3]);
incmatRR = permute(prctile(incstoRR,[2.5,50,97.5],2),[2,1,3]);

mrtmat = permute(prctile(mrtsto,[2.5,50,97.5],2),[2,1,3]);


% -------------------------------------------------------------------------
% --- Plot figure of incidence and mortality impacts ----------------------

ff=figure; lw = 1.5; fs = 14;
% allmat = cat(4,incmat,incmat2);
allmat = cat(4,incmat,incmat2,incmatRR);

cols = linspecer(size(allmat,3));
xx = [2022:2040];

lgs = {'Baseline','TPT, in-country migrants, 50% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in contacts','+ TPT in vulnerable population, 50% annually','+ ACF in migrants and vulnerable', '+ TPT in general population, 50% annually'};


plotopts = {'All incidence','RR incidence','Alternative incidence'};
plotopt  = plotopts{1};



if strcmp(plotopt, 'All incidence')
    
    % --- Single incidence plot, showing one by one
    hold on;
    for ii = 1:size(allmat,3)
        plt = allmat(:,:,ii,1);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        %jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
        
        ylim([0 7.5]);
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
        
        yline(1,'k--'); yline(0.1,'k--');
        ylabel('Incidence rate per 100,000 population');
        
        legend(lg,lgs(1:ii),'Location','Southwest');
        %pause;
    end
    
elseif strcmp(plotopt, 'RR incidence')
    
    % --- Incidence of RR-TB
    hold on;
    
    selinds = [1:3];
    allplt = squeeze(allmat(:,:,selinds,3));
    lgsel  = lgs(selinds);
    
    for ii = 1:size(allplt,3)
        plt = allplt(:,:,ii);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
        
        yline(0.1,'k--');
        ylabel('Incidence rate per 100,000 population');
        
        legend(lg,lgsel(1:ii),'Location','NorthEast');
    end
    
elseif strcmp(plotop, 'Alternative incidence')
    
    % --- Two ways of counting
    tis = {'All incidence','Incidence without treatment history'};
    for is = 1:2
        subplot(1,2,is); hold on;
    
        for ii = 1:size(allmat,3)
            plt = allmat(:,:,ii,is);
            lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
            jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
        end
        yl = ylim; yl(1) = 0; ylim(yl); ylims(is,:) = ylim;
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
    
        title(tis{is});
    
        yline(1,'k--');
        yline(0.1,'k--');
    
    end
    subplot(1,2,2);
    ylim(ylims(1,:));
    
    % legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
    legend(lg,'Baseline','TPT, in-country migrants, 50% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in vulnerable population, 50% annually','+ Find and treat active TB in migrants and vulnerable', '+ TPT in general population, 50% annually');
    subplot(1,2,1);
    ylabel('Rate per 100,000 population');

end













return;

tis = {'All incidence','Incidence without treatment history'};
for is = 1:2
    subplot(1,2,is); hold on;

    for ii = 1:size(allmat,3)
        plt = allmat(:,:,ii,is);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
    end
    yl = ylim; yl(1) = 0; ylim(yl); ylims(is,:) = ylim;
    xlim([2022 2040]);
    set(gca,'fontsize',fs);

    title(tis{is});

    yline(1,'k--');
    yline(0.1,'k--');

end
subplot(1,2,2);
ylim(ylims(1,:));

% legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
legend(lg,'Baseline','TPT, in-country migrants, 50% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in vulnerable population, 50% annually','+ ACF in migrants and vulnerable', '+ TPT in general population, 50% annually');
subplot(1,2,1);
ylabel('Rate per 100,000 population');


% --- Find remaining sources of incidence
vec = squeeze(props(end,:,end));
tmp  = abs(vec)/sum(abs(vec));
tmp2 = sortrows([tmp; 1:length(tmp)]',-1);

lbls = {}; ind = 1;
set1 = {'ds','rr'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
for is3 = 1:length(set3)
    for is2 = 1:length(set2)
        for is1 = 1:length(set1)
            lbls{ind} = [set1{is1}, ' ', set2{is2}, ' ', set3{is3}];
            ind = ind+1;
        end
    end
end

inds = tmp2(1:6,2)'
tmp2(1:6,1)'
lbls(inds)





return;

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