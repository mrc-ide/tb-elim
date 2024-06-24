function [out, aux] = get_objective(x, prm, data)

r = prm.r; p = prm.p;

r.beta  = x(1);
r.gamma = x(2);

% --- Solve the model to equilibrium
init = zeros(1,5); seed = 1e-6; init(1) = (1-seed); init(3) = seed;
geq = @(t,in)goveqs_basis2(t, in, r, p);
[t0, soln0] = ode15s(geq, [0:2e3], init, odeset('NonNegative',1:5));

% --- Find the outputs
prev = soln0(end,3)*1e5;
inc  = diff(soln0(end-1:end,5))*1e5;

% Compose the objective function
sqd = @(dat, sim) sum((1-sim./dat).^2);
out = sqd(data.prev, prev) + sqd(data.inc, inc);
out = out*100;

% Store the auxiliary outputs
aux.prev = prev; aux.inc = inc;
