function acc  = cmptAcc(resp,genMu)
mu = unique(genMu);
corrInd = genMu < 0;
r = resp;
r(corrInd) = 1 - r(corrInd);
%average across theta
for hh = 1 : numel(mu)
    ind  = find(genMu == mu(hh));
    acc(hh)        = mean(r(ind));
end
