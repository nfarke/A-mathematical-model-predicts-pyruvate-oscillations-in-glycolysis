load Resultsup4.mat

%apply forward fourier transformation and data normalization
tempnorm = Results(:,101:end) - mean(Results(:,101:end),2);
tempnorm = detrend(tempnorm,1);

[row,col] = find(tempnorm < -3);
setx = unique(col);
tempnorm(:,setx) = [];

%[row,col] = find(tempnorm > 4);
%setx = unique(col);
%tempnorm(:,setx) = [];

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

for k = 1:length(Results)
    pks = [];
    [pks,locs] = findpeaks(amplitude(k,100:end),'MinPeakheight',0.02,'MinPeakProminence',0.02);
    freq = frequency(locs);
    if ~isempty(pks)
        Out(k,1:length(freq)) = freq;
    else
        Out(k,:) = 0;
    end
end


save('FFT_info_up','Out','amplitude','frequency')



