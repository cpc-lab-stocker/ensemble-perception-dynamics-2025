function [pXgivenS_E,mapping_fun]  = pXgivenS_EC(X,S,IntN,priorMn,priorVar,Wgt,pDst,prDist) 
% X: measurments; S: stimulus in physical space
% compute p(x|s) based on efficient coding model
% IntN: internal noise in kappa 
% VarDeg: variance of the prior distribution in degree
% Wgt: weight of Gaussian component in the prior 
% exmaple call:  pXgivenS_E  = pXgivenS_EC([-pi:0.1:pi],[-pi:0.2:pi],2,16,0.99) 
snsSpc        = S; % sensory space 
stepSize      = abs(S(2) - S(1)); 
prior         = ePrior(priorMn,priorVar,Wgt,S,stepSize,pDst,prDist);  
mapping_fun   = cumtrapz(S, prior) * 2*pi-pi;                         % the cumulative mapping function
ivsStmSpc     = interp1(mapping_fun, S, snsSpc,  'linear', 'extrap'); % prior distribution an observer learned in the experiment
cnst          = 1./(2*pi*besseli(0,IntN));
pXgivenS_     = exp(IntN.*cos(X-S')).*cnst;   % p(x = 1|s) 
% asymmetric likelihood function: columns indicate likelihood functions 
pXgivenS_E    = interp1(ivsStmSpc, pXgivenS_, S, 'linear', 'extrap');
