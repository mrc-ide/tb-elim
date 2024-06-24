function [out, aux] = get_objective99(x, ref, prm, gps, calfn)

i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

[p, r] = allocate_parameters(x, p, r, xi);

% Print the parameter values and bounds
disp('Parameter values:');
disp(x);
disp('Lower bounds:');
disp(prm.bounds(1,:));
disp('Upper bounds:');
disp(prm.bounds(2,:));

tmp1  = [prm.bounds; x];
tmp2  = diff(tmp1([1,3,2],:),[],1);
cond1 = min(tmp2(:)) < 0;

if cond1
    disp('Parameter bounds violated.');
    out = -Inf;
    aux = NaN;
    return;
end

try
    M = make_model(p, r, i, s, gps);
    
    init = zeros(1, i.nx);
    seed = 1e-5;
    init(i.U.ch) = 1 - seed;
    init(i.I.ch) = seed;

    % Print initial conditions
    fprintf('Initial conditions:\n');
    disp(init);
    
    geq = @(t, in) goveqs_basis2(t, in, i, M, agg, sel, r, p);
    [t0, soln0] = ode15s(geq, [0:2020], init, odeset('NonNegative', 1:5));
    
    % Print ODE solution at initial and final steps
    fprintf('Initial ODE solution:\n');
    disp(soln0(1,:));
    fprintf('Final ODE solution:\n');
    disp(soln0(end,:));

    % Check if the ODE solution is all zeros
    if all(soln0(:) == 0)
        error('ODE solution is all zeros.');
    end

    dsol = diff(soln0, [], 1);
    incd2010 = dsol(t0 == 2010, i.aux.inc(1)) * 1e5;
    incd2020 = dsol(end, i.aux.inc(1)) * 1e5;
    incd = dsol(end, i.aux.inc) * 1e5;
    mort = dsol(end, i.aux.mort) * 1e5;
    p_adTB = incd(3) / incd(1);
    
    sfin = soln0(end, :);
    p_adpopn = sum(sfin(s.ad)) / sum(sfin(1:i.nstates));

    % Print incidences and mortality
    fprintf('Incidences and Mortality:\n');
    fprintf('incd2010: %f\n', incd2010);
    fprintf('incd2020: %f\n', incd2020);
    fprintf('incd: %f\n', incd);
    fprintf('mort: %f\n', mort);
    fprintf('p_adTB: %f\n', p_adTB);
    fprintf('p_adpopn: %f\n', p_adpopn);
    
    if any(incd < 0.01)
        disp('Incidence values are too small.');
        out = -Inf;
        aux = NaN;
    else
        out = calfn.fn(incd2010, incd2020, mort, p_adTB, p_adpopn);
        aux.soln = soln0;
        aux.incd = dsol(find(t0 == 2010):end, i.aux.inc(1)) * 1e5;
        aux.incd2010 = incd2010;
        aux.incd2020 = incd2020;
        aux.mort = mort;
        aux.p_adTB = p_adTB;
        aux.p_adpopn = p_adpopn;
    end
catch ME
    disp('Error in ODE calculation or subsequent calculations:');
    disp(ME.message);
    out = -Inf;
    aux = NaN;
end
end
