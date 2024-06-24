clear all; load calibration_res; load Model_setup;

obj  = @(x) get_objective2(x, ref, prm, gps, lhd);

ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/150);
xs  = xsto(ix0:dx:end,:);
opts = odeset('Nonnegative',1:i.nstates);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end

    xx = xs(ii,:);
    [out,aux] = obj(xx);
    init = aux.soln(end,:);

    % --- Model the baseline ----------------------------------------------
    [p0,r0] = allocate_parameters(xx, p, r, xi);
    M0 = make_model(p0,r0,i,s,gps);
    
    geq = @(t,in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
    [t,soln] = ode15s(geq, [2021:2031], init, opts);

    dsol = diff(soln,[],1);
    inc(:,ii,1)  = dsol(:,i.aux.inc(1))*1e5;
    inc1(:,ii,1) = dsol(:,i.aux.inc(4))*1e5;
    inc2(:,ii,1) = dsol(:,i.aux.inc(5))*1e5;
    inc3(:,ii,1) = dsol(:,i.aux.inc(6))*1e5;
    
    % --- Model interventions ---------------------------------------------
    p1 = p0; r1 = r0; p1.migrTPT = 1; 
    r1.TPT(1) = 0.5; 
    r1.TPT(2) = 0.5; 
    p1.TPTeff = 0.6; 
    %r1.relapse(2,3) = 0;
    r1.followup = 1e3;
    r1.ACF = 1e3;
    M1 = make_model(p1,r1,i,s,gps);
    
    geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, M1, p0, p1, [2022 2025], agg, sel, r0);
    [t,soln] = ode15s(geq, [2021:2031], init, opts);
    
    dsol = diff(soln,[],1);
    inc(:,ii,2)  = dsol(:,i.aux.inc(1))*1e5;   
    inc1(:,ii,2) = dsol(:,i.aux.inc(4))*1e5;
    inc2(:,ii,2) = dsol(:,i.aux.inc(5))*1e5;
    inc3(:,ii,2) = dsol(:,i.aux.inc(6))*1e5;
    
    % Get the break-up of incidence
    vec = dsol(end,i.aux.incsources);
    props(ii,:) = vec/sum(vec);
    
end
fprintf('\n');

incpct  = prctile(inc,[2.5,50,97.5],2);
inc1pct = prctile(inc1,[2.5,50,97.5],2);
inc2pct = prctile(inc2,[2.5,50,97.5],2);
inc3pct = prctile(inc3,[2.5,50,97.5],2);

figure; hold on;
plot(incpct(:,:,1),'b');
plot(incpct(:,:,2),'r');
plot(inc1pct(:,:,2),'g');
plot(inc2pct(:,:,2),'c');
plot(inc3pct(:,:,2),'k');
yl = ylim; yl(1) = 0; ylim(yl);

line(xlim, 0.1*[1 1], 'linestyle', '--');

props_pct = prctile(props,[2.5,50,97.5],1);
md = reshape(props_pct(2,:),2,size(props_pct,2)/2);