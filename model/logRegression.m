function loglik = logRegression(pars,data,dataType)
Beta  = pars;    % parameters of logistic regression
%lamda = para(end);
resp  = data(:,end);
x     = [ones(size(data,1),1),data(:,1:end-1)];       % independent variables (intecept,2back and 1back)
if strcmp(dataType,'linear') % Deal with circular case
    y = (1 ./ (1 + exp(-(x*Beta)))); % lamda/2+(1-lamda).*
elseif strcmp(dataType,'circular')
    w          = repmat(Beta(2:end),1,size(resp,1));
    lnCmp      = wrapToPi(circ_vecMean(data(:,1:end-1),w',2) + ones(size(data,1),1).*Beta(1));
    y          = (1 ./ (1 + exp(-lnCmp)));
end 
if numel(unique(resp)) == 2 % binary response 
    indices        = find(resp == 0);
    y(indices)     = 1 - y(indices);
    y              = max(y,1e-6); 
    loglik         = -sum(log(y));

else  % resp probability instead
    loglik =  sum(- resp.*log(y) - (1-resp).*log(1-y));
end 

end