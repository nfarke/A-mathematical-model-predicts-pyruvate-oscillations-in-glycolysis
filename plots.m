load Results10k.mat
load FFT_info_10k

PAR = PAR(:,1:22);
[p,~] = Sample(1,1);
px = p(1:22);
par_selected = PAR(find(Out(:,1)),1:22);
amp_selected = max(amplitude(find(Out(:,1)),:),[],2);

imag_selected = IMAG(find((Out(:,1))));

p={'Vmax1','Vmax2','Vmax3','Vmax4','Vmax5','Vmax6', 'k1','k2','k3','Km1','Km2',...
    'Km3','Km4','Km5','n1','n2','n3','ratio1','a1','a2','a3','a4','g6p','fbp','pep','pyr'};

%histogram plot
figure(1)
counter = 1;
for k = 7:21

nexttile
x1 = ones(size(PAR(:,k)));
x2 = ones(size(par_selected(:,k)))*2;
boxchart(x1,log10(PAR(:,k)))
hold on
boxchart(x2,log10(par_selected(:,k)))
title(p(k))
[h,px] = ttest2(log10(par_selected(:,k)),log10(PAR(:,k)),'Alpha',0.01);
P(counter) = px;
H(counter) = h;
counter = counter + 1;
end

%boxplote of the phase
figure(2)
Out(Out==0) = nan;
for k = 1:length(Out)
    freq(k,1) = max(Out(k,:));  
end

figure(2)
boxplot(1./freq,'symbol','')
hold on
r = 1.10 - 0.2*rand(length(1./freq),1);
scatter(r,1./freq,10,'filled','b')
ylabel('Period')
set(gca,'yscale','log')

idx = find(~isnan(Out(:,1)));
%plot fourier transformation / time course
% for k = 1:length(idx)
%     figure(k+2)
%     subplot(1,2,1)
%     plot(frequency,amplitude(idx(k),:))
%     subplot(1,2,2)
%     plot(1:501,Results(idx(k),:))
% end

%plot all time courses together
% for k = 1:2:10
%     figure(200)
%     plot(1:501,Results(idx(k),:))
%     hold on
% end
