clear all; 

C = readtable('TB_notifications_2022-03-11.csv');

colnames = C.Properties.VariableNames;

% --- Filter on 2019 year
rows    = find(C.year==2019);
C1      = C(rows,:);

% --- Pull out the relevant notification data
notifcols  = [find(strcmp(colnames,'iso3')), find(strcmp(colnames,'ret_nrel')), find(strcmp(colnames,'c_newinc'))];
C2         = C1(:,notifcols);
C2.allnoti = C2.ret_nrel + C2.c_newinc;

% --- Prepare for saving
notifs = C2(:,[1,end]);
save notif_data notifs;