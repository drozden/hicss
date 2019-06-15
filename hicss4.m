rng = 20:36;
nbins = 100; nt = 10000; nf = length(rng);
cov = zeros(nbins,nf); net = cov; pct = cov;
%hc = zeros(nbins,nf); hn = hc; hp = hc;
for n = 1:nf
    s = sprintf('SavedTrials%cHICSS%02d.mat',filesep,rng(n));
    d = dir(fullfile(pwd,s));
    load(fullfile(pwd,'SavedTrials',d.name));
    nt = length(outputs);
    cov(1:nt,n) = reshape([outputs.coverage],nt,1);
    net(1:nt,n) = reshape([outputs.network_traffic],nt,1);
    pct(1:nt,n) = reshape([outputs.percent_reachable],nt,1);
end
minC = min(cov(:)); maxC = max(cov(:));
minN = min(net(:)); maxN = max(net(:));
minP = min(pct(:)); maxP = max(pct(:));
binsC = linspace(minC,maxC,nbins);
binsN = linspace(minN,maxN,nbins);
binsP = linspace(minP,maxP,nbins);
hc = hist(cov,binsC);
hn = hist(net,binsN);
hp = hist(pct,binsP);
figure(6); clf; set(gcf,'Color',[1,1,1]); surf(rng,binsC,hc); 
title('Coverage'); xlabel('# Motes'); ylabel('Hist Bins'); zlabel('Freq');
figure(7); clf; set(gcf,'Color',[1,1,1]); surf(rng,binsN,hn);
title('Network Traffic'); xlabel('# Motes'); ylabel('Hist Bins'); zlabel('Freq');
figure(8); clf; set(gcf,'Color',[1,1,1]); surf(rng,binsP,hp); 
title('Percent Connected'); xlabel('# Motes'); ylabel('Hist Bins'); zlabel('Freq');
