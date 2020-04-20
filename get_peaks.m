function [exact,exact_peak]=get_peaks(r,Fs,thresh,Fmin)
	R=20*log10(r);
	s_win=(length(R)-1)*2;
	[pks,locations]= findpeaks(R,'MinPeakHeight',thresh);
	i=1;
	minpeak=Fmin;
	minpeakbin=(minpeak/Fs)*s_win+1;
	maxpeak=floor(Fs/2);
	maxpeakbin=(maxpeak/Fs)*s_win+1;
	for j=locations'
		if j<minpeakbin || j>maxpeakbin
			pks(i)=[];
			locations(i)=[];
		else
		i=i+1;
	end
	end
	exact=locations;
	exact_peak=pks;
	i=1;
	if ~isempty(locations)
	for j=locations'
		ym1=R(j-1);
		y0=R(j);
		ym2=R(j+1);
		[p,y,a]=qint(ym1,y0,ym2);
		exact(i)=j+qint(ym1,y0,ym2)-1;
		exact_peak(i)=y;
		i=i+1;
	end
	end
end
