function [Y,T] = stft(x,win,nHop,NFFT)
%% Compute the STFT
% y = STFT(x,win,nHop,win_analysis)
% params:
%   x is the input signal
%   win is the analysis window
%   nHop is the analysis hop size
%   NFFT is the number of points in each DFT
%

nWin = length(win);
L = length(x);

nFrames = floor((L-nWin)/nHop+1);
nBins = NFFT/2+1;
Y = zeros(nBins,nFrames);
T = zeros(1,nFrames);

pin = 0;
%x = [zeros(nWin/2,1);x];
x=[zeros(nWin,1);x;zeros(nWin-mod(L,nHop),1)];
for n = 1:nFrames
    grain = x(pin+1:pin+nWin).*win;
    f = fft(fftshift(grain),NFFT);
    Y(:,n) = f(1:nBins);
    T(n) = pin + 1;
    pin = pin + nHop;
end
