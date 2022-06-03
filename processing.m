load Results10k.mat

%apply forward fourier transformation and data normalization
tempnorm = Results(:,101:end) - mean(Results(:,101:end),2);
tempnorm = detrend(tempnorm,1);

[row,col] = find(tempnorm < -7);
setx = unique(row);
tempnorm(setx,:) = [];
IMAG(setx) = [];
PAR(setx,:) = [];
RE(setx) = [];
Results(setx,:) = [];
stable(setx) = [];

[row,col] = find(tempnorm > 5);
setx = unique(row);
tempnorm(setx,:) = [];
tempnorm(setx,:) = [];
IMAG(setx) = [];
PAR(setx,:) = [];
RE(setx) = [];
Results(setx,:) = [];
stable(setx) = [];

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

for k = 1:length(tempnorm)
    pks = [];
    [pks,locs] = findpeaks(amplitude(k,1000:end),'MinPeakheight',0.005,'MinPeakProminence',0.005);
    freq = frequency(locs);
    if ~isempty(pks)
        Out(k,1:length(freq)) = freq;
    else
        Out(k,:) = 0;
    end
end

save('FFT_info_10k','Out','amplitude','frequency')
save('Results10k','PAR','Results','stable','RE','IMAG','-v7.3')           



