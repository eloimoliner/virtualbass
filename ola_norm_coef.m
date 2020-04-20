function y = ola_norm_coef(win_analysis,win_synthesis,nHop_synthesis)
% Assumes constant overlap-add

nWin = length(win_analysis);
nHop = nHop_synthesis;

win = win_analysis .* win_synthesis;
idx = nWin / 2 + 1;
y = win(idx);

m = 1;
i = idx - m * nHop;
j = idx + m * nHop;
while i > 0 &&  j <= nWin
    y = y + win(i) + win(j);
    m = m + 1;
    i = idx - m * nHop;
    j = idx + m * nHop;
end
