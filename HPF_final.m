function y=HPF_final(signal,Fc,Fs)

	N=3000;  %filter order
	if (~exist('Fs', 'var'))
	Fs=44100;
    end
	if (~exist('Fc', 'var'))
	Fc=150;
    end
	SidelobeAtten=80;  % window side lobe atten
	win=chebwin(N+1,SidelobeAtten);  %chebychev window
	b=fir1(N, Fc/(Fs/2), 'high', win, 'scale');
	%figure('Name','HPF','NumberTitle','off')
	%freqz(b,1)
	HPF=dfilt.dffir(b);
	y=filter(HPF,signal);

