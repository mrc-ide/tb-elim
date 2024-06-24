% Define or load your parameters and initial conditions
i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

% Assuming allocate_parameters is a function that sets up your parameters
[p,r] = allocate_parameters(x, p, r, xi);

% Define initial conditions
init = zeros(1,i.nx);
seed = 1e-5;
init(i.U.dom)       = 1 - 0.168 - seed;
init(i.U.migr_rect) = 0.168;
init(i.I.dom)       = seed;

% Define the model and the ODE function for the initial phase
p0 = p; r0 = r; p0.betadec = 0;
r0.gamma = r0.gamma_2015;
M0 = make_model(p0, r0, i, s, gps);
geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, agg, sel, r0, p0);

% Run the ODE solver for the initial phase
[t0, soln0] = ode15s(geq0, [0:5e3], init, odeset('NonNegative',1:5));

% Define the model and the ODE function for the follow-up phase
p1 = p0; r1 = r0; r1.TPT = [0 r.TPT2020rec 0];
M1 = make_model(p1, r1, i, s, gps);
p2 = p; r2 = r; r2.gamma = r1.gamma_2020;
M2 = make_model(p2, r2, i, s, gps);
geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, [2015 2020; 2010 2020], i, s, p2, r2, prm, sel, agg);

% Run the ODE solver for the follow-up phase
[t1, soln1] = ode15s(geq1, [2010:2020], soln0(end,:), odeset('NonNegative',1:5));

% Plot the results
figure;
plot(t0, soln0);
xlabel('Time');
ylabel('State Variables');
title('Initial Simulation');

figure;
plot(t1, soln1);
xlabel('Time');
ylabel('State Variables');
title('Follow-up Simulation');
