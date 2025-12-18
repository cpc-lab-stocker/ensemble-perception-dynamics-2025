function [optPar,maxLL] = logRegFit(data,dType,mtd)
% fit logistic regression with lapse rate
nSP           = 10; 
b0            = [-1,2];
bx            = [0.1,0.5];
nItem         = size(data,2) - 1;
lb            = [-3,ones(1,nItem).* 0]; %  -.5 % constrain the weight to be non-negative 
ub            = [3,ones(1,nItem).* 10];
plb           = [-0.5,ones(1,nItem).* 0.05];
pub           = [1.5,ones(1,nItem).* 0.9];
pars          = [];    
parInd        = 1:nItem + 1;
%parInd       = 1:nItem;

%nonbcon      = []; % linear constraints on the parameters
A = []; b = [];
options = optimset('Display','iter');
options.UseParallel = true;

if strcmp(dType,'linear')
    fun  = @(x)logRegression(x,data,dType);
elseif strcmp(dType,'circ')
    fun  = @(x)logRegression(x,data,'circular');
end

for ii   = 1: nSP
    b0_ = rand(1).*(b0(2)-b0(1))+b0(1);
    bX_ = rand(nItem,1).*(bx(2)-bx(1))+bx(1);
  % bX_ = bX_./(sum(bX_));
    x0  = [b0_;bX_]; 
    if strcmp(mtd,'fmincon')
        options       = optimoptions('fmincon','Display','iter','Algorithm','sqp');
        [pars(:,ii),maxLL(ii)] = fmincon(fun,x0(parInd),A,b,[],[],lb(parInd),ub(parInd),[],options);
    elseif strcmp(mtd,'fminsearch')
        [pars(:,ii),maxLL(ii)] = fminsearch(fun,x0(parInd),options);
    elseif  strcmp(mtd,'bads')
        [pars(:,ii),maxLL(ii)] = bads(fun,x0(parInd)',lb(parInd),ub(parInd),plb(parInd),pub(parInd),[]);
    end

  %  [pars(:,ii),maxLL(ii)] = bads(fun,x0(parInd),lb(parInd),ub(parInd),plb(parInd),pub(parInd),[]);
end

[maxLL,ind]  = min(maxLL);
optPar       = pars(:,ind);

