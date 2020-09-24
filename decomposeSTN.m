function [X,Xs,Xt,Xn,Rs,Rt,Rn] = decomposeSTN(x,S,nWin,nHop,Fs)
NFFT = nWin;
%win = hann(nWin,'periodic');
win = hann(nWin);
nHopA = round(nHop/S);
%O_analysis = nWin / nHopA; % analysis overlap factor

% Compute STFT
[X,~] = stft(x,win,nHopA,NFFT);

% Compute transientness
filter_length_t = 600e-3; % in ms
filter_length_f = 180; % in Hz

nMedianH = round(filter_length_t * Fs / nHopA);
nMedianV = round(filter_length_f* NFFT / Fs);
Rt = transientness(X,nMedianH,nMedianV);

Rs = 1-Rt;
Rn = 1-sqrt(abs(Rt-Rs));
Rt=Rt-Rn/2;
Rs=Rs-Rn/2;


Xs = X.*Rs;
Xt = X.*Rt;
Xn = X.*Rn;


