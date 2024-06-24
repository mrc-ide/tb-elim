clear all;

% Natural history parameters
r.progression  = 0.0826;
r.LTBI_stabil  = 0.872;
r.reactivation = 0.0006;

r.self_cure    = 1/6;
r.relapse      = 0.003;
r.mu           = 1/66;                                                        % natural mortality
r.muTB         = 1/6;                                                         % TB related mortality
p.imm          = 0.8;                                                         % Reduced susceptibility conferred by previous infection

prm.p = p; prm.r = r; 

x = [10 .5];

r.beta  = x(1);
r.gamma = x(2);

% --- Solve the model to equilibrium
init = zeros(1,6); seed = 1e-6; init(1) = (1-seed); init(4) = seed;
geq = @(t,in)goveqs_basis2(t, in, r, p);
[t0, soln0] = ode15s(geq, [0:500], init, odeset('NonNegative',1:5));

% --- Find the outputs
prev = soln0(end,3)*1e5;
inc  = diff(soln0(end-1:end,5))*1e5;

figure; plot(t0,soln0(:,4)*1e5)



return;

data.prev = 212;
data.inc  = 257;

% --- Find the best-fitting parameters (beta, r, p) -----------------------

obj = @(x) get_objective(x, prm, data);
x1 = fminsearchbnd(obj, [10 1], [0 0], [], optimset('Display','iter'))
