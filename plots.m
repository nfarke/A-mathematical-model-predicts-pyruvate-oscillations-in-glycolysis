load Results
load FFT_info

PAR = PAR(:,1:22);
[p,~] = Sample(1,1);
px = p(1:22);
par_selected = PAR(find(Out(:,1)),1:22);

p={'Vmax1','Vmax2','Vmax3','Vmax4','Vmax5','Vmax6', 'k1','k2','k3','Km1','Km2',...
    'Km3','Km4','Km5','n1','n2','n3','ratio1','a1','a2','a3','a4','g6p','fbp','pep','pyr'};

%histogram plot
figure(1)
for k = 1:21
nexttile
histogram(log10(PAR(:,k)),200,'FaceColor',[0 0 1],'EdgeColor',[0 0 1])
hold on
histogram(log10(par_selected(:,k)),15,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])

title(p(k))
end

%boxplote of the phase
figure(2)

Out(Out==0) = nan;
for k = 1:length(Out)
    freq(k,1) = min(Out(k,:));  
end
boxplot(log(1./freq))
