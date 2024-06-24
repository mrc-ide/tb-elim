clear all; 

% -------------------------------------------------------------------------
% --- Get the estimated incidence rates for each country ------------------

C = readtable('TB_burden_countries_2022-04-12.csv');
colnames = C.Properties.VariableNames;

% --- Pull out the relevant estimate of data
cols_needed = {'country','iso3', 'year', 'e_pop_num', 'e_inc_100k'};

estimcols = zeros(1,length(cols_needed));
for ii = 1:length(estimcols)
   estimcols(ii) = find(strcmp(colnames,cols_needed{ii})); 
end
C1 = C(:,estimcols);

% --- Find which countries qualify as low-incidence
rows = find(C.year==2019 & C.e_inc_100k<=10 & C.e_pop_num>=3e5);
C1b = C(rows,:);
ctrs1 = C1b.iso3;


% -------------------------------------------------------------------------
% --- Find the contributions from the elderly -----------------------------

C = readtable('TB_burden_age_sex_2022-04-12.csv');
colnames = C.Properties.VariableNames;

% --- Pull out the relevant estimate of data
cols_needed = {'country','iso3', ...
    'sex','risk_factor','age_group','best'};

estimcols = zeros(1,length(cols_needed));
for ii = 1:length(estimcols)
   estimcols(ii) = find(strcmp(colnames,cols_needed{ii})); 
end
C2 = C(:,estimcols);


% --- For each country, find the proportion of incidence coming from >65s

ctrs2 = unique(C2.iso3);

for ic = 1:length(ctrs2)
    
    ctr = ctrs2{ic};
   if ismember(ctr, ctrs1) 
   rows = find(strcmp(C2.iso3,ctrs2{ic})); 
   Ct = C2(rows,:);
   
   % Get the total count
   row = find(strcmp(Ct.sex, 'a') & strcmp(Ct.risk_factor, 'all') & strcmp(Ct.age_group, 'all'));
   den = Ct.best(row);
   
   % Get the numbers for >65s
   row = find(strcmp(Ct.risk_factor, 'all') & strcmp(Ct.age_group, '65plus'));
   num = sum(Ct.best(row));
   
   % Calculate the ratio
   rat(ic) = num/den;
   
   else
      rat(ic) = nan; 
   end
end
% C2.rat = rat;
rat = rat';

tbl  = table(ctrs2,rat);

mat  = sortrows([rat'; 1:length(rat)]',-1);
take = find(~isnan(mat(:,1)));
ord  = mat(take,:);

fprintf('Countries with the highest old-age contributions:\n');
disp(tbl(ord(1:10,2),:));



% -------------------------------------------------------------------------
% --- Find the contributions from the foreign-born

C = readtable('TB_notifications_2022-04-12.csv');
colnames = C.Properties.VariableNames;

% --- Pull out the relevant estimate of data
cols_needed = {'country','iso3', 'year', 'notif_foreign', 'ret_nrel', 'c_newinc'};

estimcols = zeros(1,length(cols_needed));
for ii = 1:length(estimcols)
   estimcols(ii) = find(strcmp(colnames,cols_needed{ii})); 
end
C3  = C(:,estimcols);
C3b = C3(C3.year==2019,:);
C3b.prop = C3b{:,4}./sum(C3b{:,[5,6]},2);

% Filter on the low-incidence countries
rows = find(ismember(C3b.iso3, ctrs1));
C3c = C3b(rows,:);

ctrs3 = C3c.iso3;

props = C3c.prop;
mat   = sortrows([props'; 1:length(props)]',-1);
take  = find(~isnan(mat(:,1)));
ord   = mat(take,:);

fprintf('Countries with the highest foreign contributions: \n');
rows = ord(1:10,2);
disp(C3c(rows,[2,end]))

