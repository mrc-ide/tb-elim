function out = goveqs_scaleup(t, in, i, M0, M1, times, agg, sel, r, p)

scale = max(min((t-times(1))/(times(2)-times(1)),1),0);
Ms = M1;
Ms.lin = M0.lin + scale*(M1.lin-M0.lin);

out = goveqs_basis2(t, in, i, Ms, agg, sel, r, p);