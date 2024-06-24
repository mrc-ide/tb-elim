load analysis_outputs.mat

figure; histogram(ch_inc)
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (per 100,000 children-years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');


[yy, xx] = ecdf(ch_inc);
%[yy2, xx2] = ecdf(ch_inc2);

figure('Position', [577, 226, 1029, 732]); hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 
xlim([min(xx), 0.18]);
ylim([0, 1]);
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');

boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));
boundary_90 = xx(find(yy >= 0.1, 1));
boundary_40 = xx(find(yy >= 0.6, 1));

line([boundary_90, boundary_90], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([boundary_40, boundary_40], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

fill([min(xx); xx(xx <= boundary_90); boundary_90], [0; 1-yy(xx <= boundary_90); 0], [1, 0.84, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Gold'); % Gold
fill([boundary_90; xx(xx > boundary_90 & xx <= boundary_40); boundary_40], [0; 1-yy(xx > boundary_90 & xx <= boundary_40); 0], [0.75, 0.75, 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Silver'); % Silver
fill([boundary_40; xx(xx > boundary_40); max(xx)], [0; 1-yy(xx > boundary_40); 0], [0.8, 0.5, 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Bronze'); % Bronze

% Add horizontal lines at y = 0.90 and y = 0.40
line([0, boundary_90], [1-0.1, 1-0.1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([0, boundary_40], [1-0.6, 1-0.6], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% Add filled circle dots where the lines meet
plot(boundary_90, 1-0.1, 'ko', 'MarkerFaceColor', 'k');
plot(boundary_40, 1-0.6, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); 

xlim([min(xx), 0.18]);
%title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');
legend('Location', 'best', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');
hold off;






%%%%%%

% get data 
[yy, xx] = ecdf(ch_inc2);
%[yy2, xx2] = ecdf(ch_inc2);


figure('Position', [577, 226, 1029, 732]); hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 

% get bounds
boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));

% plot verticals
line([boundary_66, boundary_66], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([boundary_33, boundary_33], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% shade
fill([min(xx); xx(xx <= boundary_33); boundary_33], [0; 1-yy(xx <= boundary_33); 0], [1, 0.84, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Gold'); % Gold
fill([boundary_33; xx(xx > boundary_33 & xx <= boundary_66); boundary_66], [0; 1-yy(xx > boundary_33 & xx <= boundary_66); 0], [0.75, 0.75, 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Silver'); % Silver
fill([boundary_66; xx(xx > boundary_66); max(xx)], [0; 1-yy(xx > boundary_66); 0], [0.8, 0.5, 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Bronze'); % Bronze

%  x-axis limit
%xlim([min(xx), 0.18]);

title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (cases per 100,000 children)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');

legend('Location', 'best', 'FontWeight', 'bold');

set(gca, 'FontWeight', 'bold');

hold off;



[yy1, xx1] = ecdf(ch_inc);
[yy2, xx2] = ecdf(ch_inc2);


figure('Position', [577, 226, 1029, 732]); hold on;
plot(xx1, 1-yy1, 'Color', [0/255, 0/255, 0/255], 'LineWidth', 2); %  blue
plot(xx2, 1-yy2, 'Color', [255/255, 0, 0], 'LineWidth', 2); %  red


ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (cases per 100,000 children)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontSize', 14, 'FontWeight', 'bold');
title('Thresholds for reaching elimination vs pre-elimination', 'FontSize', 16, 'FontWeight', 'bold');

hold off;


%tornado plot


vec = pcm(end, 1:end-1);


param_names = {'Contact rate', 'Secular trend in contact rate', 'Per-capita rate of treatment initiation', 'Child mortality', 'Relapse risk', 'Ageing', 'Relative infectivity, children vs adults', 'Intergenerational mixing', 'Baseline incidence'};


[~, sorted_indices] = sort(abs(vec), 'descend');
sorted_vec = vec(sorted_indices);


sorted_param_names = param_names(sorted_indices);

sorted_vec = flipud(sorted_vec);
sorted_param_names = flipud(sorted_param_names);


figure;
barh(sorted_vec, 'red');


set(gca, 'YDir', 'reverse');


xlabel('Partial Correlation', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Parameters', 'FontSize', 14, 'FontWeight', 'bold');
title('Partial Correlation Ordered by Absolute Value', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'YTick', 1:length(sorted_param_names), 'YTickLabel', sorted_param_names);


set(gca, 'FontWeight', 'bold');




%%%%%%

% re do no shading with bounds




[yy, xx] = ecdf(ch_inc);
%[yy2, xx2] = ecdf(ch_inc2);

figure('Position', [577,   226 ,   1029 ,732]); hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 
xlim([min(xx), 0.2]);
ylim([0, 1]);
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
%title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (cases per 100,000 children)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');

boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));
boundary_90 = xx(find(yy >= 0.1, 1));
boundary_40 = xx(find(yy >= 0.6, 1));

line([boundary_90, boundary_90], [0, 1-0.1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% Add horizontal lines at y = 0.66 and y = 0.33
line([0, boundary_90], [1-0.1, 1-0.1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% Add filled circle dots where the lines meet
plot(boundary_90, 1-0.1, 'ko', 'MarkerFaceColor', 'k');

xlim([min(xx), 0.18]);
%title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');
%legend('Location', 'best', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');
hold off;


%%elim pre elim shaded
[yy1, xx1] = ecdf(ch_inc);
[yy2, xx2] = ecdf(ch_inc2);

figure('Position', [577, 226, 1029, 732]); hold on;

% Plot the ECDF curves
plot(xx1, 1-yy1, 'Color', [0/255, 0/255, 0/255], 'LineWidth', 2); % blue
plot(xx2, 1-yy2, 'Color', [255/255, 0, 0], 'LineWidth', 2); % red

% Fill the Bronze region (between blue and red curves)
fill([xx2; flip(xx1)], [1-yy2; flip(1-yy1)], [205/255, 127/255, 50/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Define x-axis regions for silver and gold
silver_x = linspace(0.05, 0.2, 100); % Region from 0.05 to 0.1
gold_x = linspace(0, 0.05, 100); % Region from 0 to 0.05

% Ensure unique x-values for interpolation
[xx1_unique, idx1] = unique(xx1);

% Interpolate to find y-values for the blue curve
silver_y1 = interp1(xx1_unique, 1-yy1(idx1), silver_x, 'linear', 'extrap');
gold_y1 = interp1(xx1_unique, 1-yy1(idx1), gold_x, 'linear', 'extrap');

% Ensure that the lengths match for fill function
if length(silver_x) == length(silver_y1)
    % Fill the Silver region
    fill([silver_x, flip(silver_x)], [zeros(size(silver_x)), flip(silver_y1)], [192/255, 192/255, 192/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
else
    disp('Lengths of x and y vectors for the silver region do not match.');
end

if length(gold_x) == length(gold_y1)
    % Fill the Gold region
    fill([gold_x, flip(gold_x)], [zeros(size(gold_x)), flip(gold_y1)], [255/255, 215/255, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
else
    disp('Lengths of x and y vectors for the gold region do not match.');
end

% Set axis properties
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
ylim([0,1])
xlabel('Incidence in children (cases per 100,000 children-years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability that country has reached (pre) elimination', 'FontSize', 14, 'FontWeight', 'bold');
title('Thresholds for reaching elimination vs pre-elimination', 'FontSize', 16, 'FontWeight', 'bold');

hold off;

