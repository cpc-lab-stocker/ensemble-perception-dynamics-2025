function respProb = eftMdl(IntN,Wgt,data,dsetMu,setSigma,pDst,prDist)
sc          = 2;
stepSize    = deg2rad(1);
stmSpc      = -pi:stepSize:pi;
X           = stmSpc; % measurements
nMsmSmps    = 2000; 
setSize     = size(data,1);
nTr         = size(data,2);
dataRad     = deg2rad(data).* sc; 
priorVar    = setSigma.* sc;
priorMn     = dsetMu.* sc; 
[pX_igivenS_i,mapping_fun]  = pXgivenS_EC(X,stmSpc,IntN,priorMn,priorVar,Wgt,pDst,prDist);
msmSmps     = nan(nMsmSmps,setSize,nTr);
for jj = 1 : nTr
    for ii  = 1:setSize
        STilde               = interp1(stmSpc, mapping_fun, dataRad(ii,jj), 'linear', 'extrap');
        msmSmps(:,ii,jj)     = vmrand(STilde, IntN, [nMsmSmps,1]);
    end
end
%
% fixed across trials: likelihood functions
S_i                   = X;
mu_j1  = deg2rad(dsetMu(:,1).*sc); mu_j0 = deg2rad(dsetMu(:,2).*sc);
sigmaC = deg2rad(setSigma.*sc); sigma_k1 = sigmaC; sigma_k0 = sigmaC;
%pGaussPar_jkgivenC1   = 1/(numel(mu_j1) * numel(sigma_k1));  
%pGaussPar_jkgivenC0   = 1/(numel(mu_j0) * numel(sigma_k0));
pGaussPar_jkgivenC1   = pDst(1:numel(mu_j1))./numel(sigma_k1);  
pGaussPar_jkgivenC0   = pDst(numel(mu_j0)+1:end)./numel(sigma_k0);

pMPar1                = nan(numel(S_i),numel(sigma_k1),numel(mu_j1));
pMPar0                = nan(numel(S_i),numel(sigma_k1),numel(mu_j1));
for ii = 1: numel(mu_j1)
    for jj = 1 : numel(sigma_k1)
        pS_igivenGaussPar_jk1               = normpdf(S_i, mu_j1(ii),sigma_k1(jj));
        pS_igivenGaussPar_jk0               = normpdf(S_i, mu_j0(ii),sigma_k0(jj)); 
        pMPar1(:,jj,ii) = pX_igivenS_i'*pS_igivenGaussPar_jk1'.*pGaussPar_jkgivenC1(ii);
        pMPar0(:,jj,ii) = pX_igivenS_i'*pS_igivenGaussPar_jk0'.*pGaussPar_jkgivenC0(ii);
    end
end

respProb = nan(nTr,1);
for gg = 1: nTr    
    pXgivenGaussPar1      = nan(numel(sigma_k1),numel(mu_j1),nMsmSmps);
    pXgivenGaussPar0      = nan(numel(sigma_k1),numel(mu_j1),nMsmSmps);  
    for ii = 1: numel(mu_j1)
        for jj = 1 : numel(sigma_k1) 
            msms1     = interp1(S_i, pMPar1(:,jj,ii), msmSmps(:,:,gg),  'linear', 'extrap'); 
            msms0     = interp1(S_i, pMPar0(:,jj,ii), msmSmps(:,:,gg),  'linear', 'extrap'); 
            pXgivenGaussPar1(jj,ii,:) = prod(msms1,2);
            pXgivenGaussPar0(jj,ii,:) = prod(msms0,2);
        end 
    end      

    pXgivenC1 = sum(pXgivenGaussPar1,[1,2]);
    pXgivenC0 = sum(pXgivenGaussPar0,[1,2]);
    
    pC1 = 0.5; pC0 = 0.5;
    pC1givenX = pXgivenC1.*pC1;
    PC0givenX = pXgivenC0.*pC0;
    d         = pC1givenX./PC0givenX; % decision variable, posterior ratio
    d         = d > 1;
    respProb(gg)  = mean(d);
end

end 