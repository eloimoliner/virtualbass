function y = istft(X,nHop,win,win_analysis)
%% Compute the inverse STFT
% y = ISTFT(X,nHop,win,win_analysis)
% params:
%   X is input spectrogram (positive frequencies and DC)
%   nHop is the hop size
%   win is the synthesis window (default is rectangular)
%   win_analysis is the analyis window used for computing spectrogram X
%

nWin = length(win);
[nBins,nFrames] = size(X);
NFFT = (nBins-1)*2;

% Window function
if nargin < 3
    % Default synthesis window
    win = rectwin(NFFT);
end
nWin = length(win);
if nargin < 4
    % Default analysis window
    win_analysis = rectwin(nWin);
end

% Length of output
L = (nFrames-1)*nHop + nWin;
y = zeros(L,1);

% OLA normalization
norm_coef = ola_norm_coef(win_analysis,win,nHop);

% Compute two-sided spectrogram
XF = zeros(NFFT,nFrames);
XF(1:nBins,:) = X(1:nBins,:);
XF(nBins+1:end,:) = conj(flipud(X(2:end-1,:)));

% Overlap-add synthesis
p = 0;
for n = 1:nFrames
    grain = fftshift(real(ifft(XF(:,n))));
    grain = grain(1:nWin).*win ./ norm_coef;
    y(p+1:p+nWin) = y(p+1:p+nWin) + grain;
    p = p + nHop;
end
