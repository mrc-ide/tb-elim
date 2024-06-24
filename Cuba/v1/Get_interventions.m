clear all; load calibration_res5.mat; load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, lhd);
ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/200);
xs  = xsto(ix0:dx:end,:);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);



timereached = NaN(size(xs, 1), 1); % empty for time reached
ch_inc = NaN(size(xs, 1), 1); % 


mk = round(size(xs, 1) / 25);
for ii = 1:size(xs, 1)
    
    if mod(ii, mk) == 0
        fprintf('%0.5g ', ii / mk); 
    end
    
    xx = xs(ii, :);
    [out, aux] = obj(xx);
    
    init = aux.soln(end, :);

    [p0, r0] = allocate_parameters(xx, p, r, xi);
    M0 = make_model(p0, r0, i, s, gps);

    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 3.91202 * [1 1];
    r4.ACF = 3.91202 * [1 1];
    M4 = make_model(p4, r4, i, s, gps); % ACF and TPT in everyone
    
    models = {M0, M4};
    
    for mi = 1:length(models)
        geq = @(t, in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t, soln] = ode15s(geq, [2022:2100], init, opts);
        
        sdiff = diff(soln, [], 1);
        
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1));
        incsto2(:,ii,mi) = sdiff(:,i.aux.inc(2));
        incsto3(:,ii,mi) = sdiff(:,i.aux.inc(3));

        % find  year when inc reaches 1 per mill
        idx = find(sdiff(:, i.aux.inc(1)) * 1e5 <= 0.1, 1); % 1 per million = 0.1 per 100,000
        if ~isempty(idx)
            timereached(ii, mi) = t(idx);
            if mi == 2 %  model 4 is second model in the list
                ch_inc(ii, 1) = soln(idx, i.aux.inc(2)); % *1e5 for per 100,000 children
            end

        else
            timereached(ii, mi) = NaN; % Set to NaN if threshold is never reached
        end
    end
end

fprintf('\n');

samplescell = mat2cell(xs, ones(size(xs, 1), 1), size(xs, 2));
resultstable = table(samplescell, timereached, ch_inc, 'VariableNames', {'Parameter Set', 'Year reached', 'Child incidence at this point in time'});
disp(resultstable);

incsto = incsto*1e5;
incsto2 = incsto2*1e5;
incsto3 = incsto3*1e5;

%parameter set that gives latest elimination
figure;
subplot(1,3,1);
plot(2022:2099, squeeze(incsto2(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); %children

subplot(1,3,2);
plot(2022:2099, squeeze(incsto3(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); %adults

subplot(1,3,3);
plot(2022:2099, squeeze(incsto(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); % all incidence



%  1st and last year of elimination
[earliest_yr, earliest_idx] = min(timereached(:, 2));
[latest_yr, latest_idx] = max(timereached(:, 2));
% params
earliest_params = xs(earliest_idx, :);
latest_params = xs(latest_idx, :);

disp('Earliest Year Parameter Set:');
disp(earliest_params);
disp(['Earliest Year: ', num2str(earliest_yr)]);

disp('Latest Year Parameter Set:');
disp(latest_params);
disp(['Latest Year: ', num2str(latest_yr)]);

% compare
comp_tab = table(earliest_params', latest_params', 'VariableNames', {'Earliest Year Params', 'Latest Year Params'});
disp('Comparison of Parameter Sets:');
disp(comp_tab);

return;



%


















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
    M2 = make_model(p2,r2,i,s,gps); %acf in every1

    p3 = p0; r3 = r0;
    p3.migrTPT = 1;
    r3.ACF = 0.69*[1 1];
    r3.TPT = 0.69*[0 1];
    M3 = make_model(p3,r3,i,s,gps); % tpt in adults
    
    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 0.69*[1 1];
    r4.ACF = 0.69*[1 1];
    M4 = make_model(p4,r4,i,s,gps); % acf and tpt in every1
    
    models = {M0, M4};
    
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t,soln] = ode15s(geq, [2022:2041], init, opts);
        
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

xx = [2022:2040];
for ii = 1:size(mat,3)
   plt = mat(:,:,ii);
   lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
   jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
end
yl = ylim; yl(1) = 0; ylim(yl);
set(gca,'fontsize',fs);
ylabel('Incidence per 100,000 population');
xlabel('Year');
yline(0.1,'k--','LineWidth', 2);
yline(1,'k--','LineWidth', 2);
title('Sampled params');

legend(lg, 'Baseline','ACF + TPT in everyone','location','SouthWest');


return;
%




paramset = xs(88, :);



[out, aux] = obj(paramset);
init = aux.soln(end, :);

[p0, r0] = allocate_parameters(paramset, p, r, xi);
M0 = make_model(p0, r0, i, s, gps);

% ---------------------------------------------------------------------
% --- Model baseline
% geq = @(t, in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
% [t, soln] = ode15s(geq, [2022:2031], init);
% sdiff = diff(soln, [], 1);
% incsto(:, 1, 1) = sdiff(:, i.aux.inc(1)) * 1e5;

% ---------------------------------------------------------------------
% --- Model intervention

p1 = p0; r1 = r0;
p1.migrTPT = 1;
M1 = make_model(p1, r1, i, s, gps);

p2 = p0; r2 = r0;
p2.migrTPT = 1;
r2.ACF = 0.69 * [1 1];
M2 = make_model(p2, r2, i, s, gps); % ACF in everyone

p3 = p0; r3 = r0;
p3.migrTPT = 1;
r3.ACF = 0.69 * [1 1];
r3.TPT = 0.69 * [0 1];
M3 = make_model(p3, r3, i, s, gps); % TPT in adults

p4 = p0; r4 = r0;
p4.migrTPT = 1;
r4.TPT = 0.69 * [1 1];
r4.ACF = 0.69 * [1 1];
M4 = make_model(p4, r4, i, s, gps); % ACF and TPT in everyone

models = {M0, M4};

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);

%incsto = zeros(2200 - 2022, 1, length(models)); 
%props = zeros(1, length(i.aux.incsources), length(models));

for mi = 1:length(models)
    geq = @(t, in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
    [t, soln] = ode15s(geq, [2022:2100], init, opts);

    sdiff = diff(soln, [], 1);
    incsto(:, 1, mi) = sdiff(:, i.aux.inc(1)) * 1e5;


    vec = sdiff(end, i.aux.incsources) * 1e5;
    props(1, :, mi) = vec / sum(vec);
end

cols = linspecer(length(models));
figure; lw = 3; fs = 14;

xx = 2022:2100-1;
for ii = 1:length(models)
    plt = incsto(:, 1, ii);
    lg(ii, :) = plot(xx, plt, 'Color', cols(ii, :), 'linewidth', lw); hold on;
end

yl = ylim; yl(1) = 0; ylim(yl);
set(gca, 'fontsize', fs);
ylabel('Incidence per 100,000 population');
xlabel('Year');
yline(0.1, 'k--', 'LineWidth', 2);
yline(1, 'k--', 'LineWidth', 2);

legend(lg, 'Baseline', 'ACF in everyone', 'ACF + adults TPT', 'ACF + TPT in everyone', 'location', 'SouthWest');
title('Single parameter set');




%




