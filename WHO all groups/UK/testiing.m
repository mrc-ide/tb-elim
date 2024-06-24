clear all; load ../v0/maxes;

load Model_setup; 

prm.p.RR_inmigr = 0;
prm.r.RR_acqu   = 0;

x1 = [xmax(1), 0, xmax(2:end)];

obj  = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

[out, aux] = obj(x1);

inds = setdiff(1:i.nstates,[s.rr,s.Tx2]);

% tmp = aux.M0.nlin.dom.ds + aux.M0.nlin.migr_rect.ds + aux.M0.nlin.migr_long.ds;
% Mnew = tmp(inds,inds);
% Mold = auxmax.M0.nlin;

tmp  = aux.M0.lam(1,:);
Mnew = tmp(:,inds);
Mold = auxmax.M0.lam;


Mnew = aux.M0.mort(inds,:);
Mold = auxmax.M0.mort;



[row, col] = find(Mnew~=Mold)
