% Efficient ensemble coding model 
%% data simulation
clear; %close all
% simulate expt data
dsetMu        = [5 10 20;-5 -10 -20]';       % generative model mean relative to the reference
setSigma      = 10; %[10,10];                % variance
p             = ones(1,numel(dsetMu(:)))/(numel(dsetMu(:)));       %  Gaussian
dDist         = 'gauss';    
[data,genMu]  = smpData(dsetMu,setSigma,12,p,1000,0);
%% model prediction
% parameters
kappa = 1;   % sensory noise 
beta  = 0.6; % strength of efficient coding
pcw   = eftMdl(kappa,beta,data,dsetMu,setSigma,p,dDist); % predict choice probability
%% analysis
% number of bins 
nBin                  = 4;
[binX,edges]          = tallyIntoAbsBins(data,nBin,40,[],'linear');
% predicted accuracy 
acc                   = cmptAcc(pcw,genMu);
% extracted weighting kernel
binCoefs              = logRegFit([binX',pcw],'linear','fmincon');
%% plotting
accClp = (acc([3,2,1]) + acc([4,5,6]))./2; 
% plot accuracy 
figure
tiledlayout(1,1,"Padding","tight","TileSpacing","tight")
set(gcf,'Position',[500,500,320,300]);
plot(abs(dsetMu),accClp','k-','LineWidth',3);
hold on % 
formatFigure('|\mu| (deg)','Accuracy')
xticks([5 10 20])
ylim([0.5,1])
%set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.04])
box off

% plot weighting kernal 
figure
tiledlayout(1,1,"Padding","tight","TileSpacing","tight")
set(gcf,'Position',[500,500,320,300]);
plot(5:10:35,binCoefs(2:end),'k-','LineWidth',3);
hold on % 
formatFigure('Relative orientation (deg)','Weight')
xticks(5:10:35)
ylim([0,2])
%set(gca,'TickDir','out');
set(gca,'TickLength',[0.02, 0.04])
box off
