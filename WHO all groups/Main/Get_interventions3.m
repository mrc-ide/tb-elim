clear all; load calibration_res_isa.mat; 

obj = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

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
    
    TPTcov = -log(1e-8);
    TPTcov = 100;

    % % Better treatment outcomes overall
    % ra = r0; pa = p0;
    % %ra.relapse(1) = 0;
    % Ma = make_model(pa,ra,i,s,gps,prm.contmat);

    % Elimination in migrants
    rb = r0; pb = p0;
    pb.TPTeff = 1;
    rb.TPT = TPTcov*[0 1 1 0];
    rb.migrTPT = 1;
    cont2 = prm.contmat; %cont2(2,:) = 0;
    Mb = make_model(pb,rb,i,s,gps,cont2);

    % Elimination in vulnerables
    rc = rb; pc = pb;
    rc.TPT = TPTcov*[0 1 1 1];
    cont3 = cont2; %cont3(3,:) = 0;
    Mc = make_model(pc,rc,i,s,gps,cont3);
    
    % Elimination in DR-TB
    rd = rc; pd = pc;
    rd.RR_acqu = 0; pd.relbeta_RR = 0;
    Md = make_model(pd,rd,i,s,gps,cont3);

    % Elimination in domestic
    re = rd; pe = pd;
    re.TPT = TPTcov*[1 1 1 1];
    Me = make_model(pe,re,i,s,gps,cont3);

    % Elimination of recurrence
    rf = re; pf = pe;
    rf.relapse(1:3) = 0;
    Mf = make_model(pf,rf,i,s,gps,cont3);

    models = {M0, Mb, Mc, Md, Me, Mf};    
    
    for mi = 1:length(models)
        
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pa, [2024 2029], agg, sel, r0);
        [t,soln] = ode15s(geq, [2022:2041], init);
        
        endsolsto(mi,:) = soln(end,:);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        mrtsto(:,ii,mi) = sdiff(:,i.aux.mort)*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;
        props(ii,:,mi) = vec/sum(vec);
        vec1 = vec;
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
xx = [2022:2040];

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
yline(1,'k--');
yline(0.1,'k--');
% legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
legend(lg,'0','a','b','c','d');
ylabel('Rate per 100,000 population');


% --- Find remaining sources of incidence
vec = squeeze(props(end,:,end));


                                              % 1.DS/RR, 2.Dom/Migr/Vuln, 3.L/P/R/T


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


