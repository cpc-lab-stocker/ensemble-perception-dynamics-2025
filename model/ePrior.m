function stmPrior = ePrior(mu,theta,w,stmSpc,stepSize,pDst,prDist)
% specify the prior distribution, which is a mixture of
% Gaussian and uniform distributions
% mu    = [5,10,15,20 25;-5,-10,-15,-20 -25]';       % set mean relative to the reference
% mu    = [10,15,20 25;-10,-15,-20 -25]';       % set mean relative to the reference

% theta = [16,16]./2;  
% p     = [0.08,0.10,0.12,0.20,0.08,0.10,0.12,0.20];  % bimodal                       % probability of sampling from diff distributions
% p     = [0.15,0.15,0.08,0.06,0.06,0.15,0.15,0.08,0.06,0.06];  % gauss
% example call: stmPrior = ePrior(mu,theta,1,[-pi:0.1:pi],0.1,p,'biMod',1);

if ~exist('bPLOT','var')       || isempty(bPLOT)       bPLOT = 0;             end
if ~exist('stepSize','var')    || isempty(stepSize)    stepSize = abs(stmSpc(2)-stmSpc(1)); end
if ~exist('mu','var')          || isempty(mu)          mu = 0; end
if ~exist('theta','var')       || isempty(theta)       theta = 0; end
if ~exist('pDst','var')        || isempty(pDst)        pDst = 1/(numel(mu(:)) * numel(theta(:))); end

unifCpt        = 1./range(stmSpc);
muRad          = deg2rad(mu);
thetaRad       = deg2rad(theta); 
if strcmp(prDist,'gauss')    
    muRad_     = repmat(muRad(:),1,numel(thetaRad(:)));
    sigmaRad_  = repmat(thetaRad,numel(muRad(:)),1); 
    pDst_      = repmat(pDst',1,numel(thetaRad(:)));
    pDst_      = pDst_./sum(pDst_,'all'); 
    mixGauss   = sum(pDst_(:).*normpdf(stmSpc,muRad_(:),sigmaRad_(:)));
end 
stmPrior       = mixGauss.*w + unifCpt.*(1-w);  % uniform distribution 