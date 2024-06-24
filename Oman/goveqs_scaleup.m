function out = goveqs_scaleup(t, in, i, s, M0, M1, p0, p1, times, agg, sel, r)

scale = min(max((t-times(1))/(times(2)-times(1)),0),1);

Ms = M1; 

Ms.lin = M0.lin     + scale*(M1.lin-M0.lin);
Ms.migrentries = M0.migrentries + scale*(M1.migrentries-M0.migrentries);

ps = p1; ps.migrTPT = p0.migrTPT + scale*(p1.migrTPT-p0.migrTPT);          % Specified separately since this enters through governing equations, not model matrix

out = goveqs_basis3(t, in, i, s, Ms, agg, sel, r, ps);