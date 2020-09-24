function w=hann(M)
w = .5*(1 - cos(2*pi*(0:M-1)'/(M-1)));
