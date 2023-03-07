rng(6)
N = 3;
config = dec2bin(0:2^N-1)' - '0';
config = config';
%parpool(17)
for k = 1:length(config)
     M = create_SS_solutions(config(k,:),num2str(k));
end
for k = 1:length(config)
    results = load(strcat('Results',num2str(k),'.mat'));
    Results = results.Results;
    subplot(1,8,k)
    plot(Results')
end

for k = 1:length(config)
    results = load(strcat('Results',num2str(k),'.mat'));
    Results = results.Results;
    tempnorm = Results(:,50:end) - mean(Results(:,50:end),2);
    %tempnorm = detrend(tempnorm,2);
    plot(tempnorm')
    
    Fs = 1;                    % Sampling frequency
    T = 1/Fs;                  % Sampling period
    L = length(tempnorm);      % Length of signal
    t = (0:L-1)*T;  

    n = 2^nextpow2(length(tempnorm));
    dim = 2;
    Y = fft(tempnorm,n,dim);
    P2 = abs(Y/L);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);

    frequency = 0:(Fs/n):(Fs/2-Fs/n);
    amplitude = P1(:,1:n/2);

    for j = 1:size(tempnorm,1)
        pks = [];
        [pks,locs] = findpeaks(amplitude(j,50:end),'MinPeakheight',0.001,'MinPeakProminence',0.001);
        freq = frequency(locs);
        if ~isempty(pks)
            Out(j,1:length(freq)) = freq;
        else
            Out(j,:) = 0;
        end
    end
    idx = find((Out(:,1)));

    
    COUNT(k,1) = length(idx);
end
