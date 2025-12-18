function smps = sampler(mu,sigma,nSmp,smpMtd)
smps          = [];
for ii = 1 : nSmp(2)
    smps(:,ii) = normrnd(mu,sigma,nSmp(1),1);
    if strcmp(smpMtd,'genSts')
        while abs(mean(smps(:,ii) - mu)) > 0.5 || abs(std(smps(:,ii)) - sigma) > 2
            smps(:,ii) = normrnd(mu,sigma,nSmp(1),1);
        end
    elseif strcmp(smpMtd,'genMu')
        while abs(mean(smps(:,ii) - mu)) > 0.5 
            smps(:,ii) = normrnd(mu,sigma,nSmp(1),1);
        end       
    end
end
