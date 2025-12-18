function [vBins,edges] = tallyIntoAbsBins(X,nBin,sp,edges,dType)
if  ~exist('edges','var') || isempty(edges)
    if strcmp(dType,'circ')
        sp = sp * 2;
    end 
    edgesPos         = linspace(0,sp,nBin + 1)./180.*pi;
    edgesPos(1)      = 0; edgesPos(end) = max(X,[],'all')./180.*pi;
    edgesNeg         = linspace(-sp,0,nBin + 1)./180.*pi;
    edgesNeg(1)      = min(X,[],'all')./180.*pi; edgesNeg(end) = 0;
end
nTr               = length(X);
X                 = X./180.*pi;
bin_Pos           = discretize(X,edgesPos,'IncludedEdge','left');
bin_Neg           = 5 - discretize(X,edgesNeg,'IncludedEdge','left');
bin_Pos(isnan(bin_Pos)) = 0; bin_Neg(isnan(bin_Neg)) = 0; 
bin_              = bin_Pos + bin_Neg; 
vBins             = nan(nBin,nTr);

for jj = 1 : nTr
    for ii = 1:nBin
        idx    =  find(bin_(:,jj) == ii);
        if ~isempty(idx)
            if strcmp(dType,'linear')
                vBins(ii,jj) = sum(X(idx,jj));
            elseif strcmp(dType,'circ')
                vBins(ii,jj) = wrapToPi(sum(X(idx,jj)));
                
            end
        else
            vBins(ii,jj) = 0; % replace with the mean of the range
        end
    end
end





