function [smp,genMu,genTheta,smpMu,smpTheta] = smpData(mu,theta,nSmp,pDist,nTrl,bPlt)
smp = []; genMu = []; genTheta = [];
smpMtd   = [];         % 'genSts''smpSts'; % sample vs. generative statistics (mean and std)
for jj   = 1 : numel(mu(:))
    for ii = 1 : numel(theta(:))
        smp_        = sampler(mu(jj),theta(ii),[nSmp,round(nTrl.*pDist(jj))],smpMtd);
        smp         = [smp_,smp];
        genMu       = [repmat(mu(jj),size(smp_,2),1);genMu];
        genTheta    = [repmat(theta(ii),size(smp_,2),1);genTheta];
    end
end
% sample mean and std
smpMu    = mean(smp); 
smpTheta = std(smp);

if bPlt == 1
    figure('position',[680   666   1000   280]);
    subplot(1,3,1)
    histogram(smp(:))
    subplot(1,3,2)
    histogram(genMu)
    formatFigure([],[]);
    subplot(1,3,3)
    histogram(genTheta)
    formatFigure([],[]);
end