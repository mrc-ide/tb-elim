clear all; load calibration_res;

ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/150);
xs  = xsto(ix0:dx:end,:);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    [out,aux] = obj(xs(ii,:));
    dsol = diff(aux.soln,[],1);
    inc(:,ii) = dsol(:,i.aux.inc(1))*1e5;
end
fprintf('\n');

incpct = prctile(inc,[2.5,50,97.5],2)*1e5;
figure; plot(incpct)