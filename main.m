clear all
close all
%%%Load audio sample
%
Fs=44100;
%audio=audioread('audios/jazz.wav');
%%audio=audioread('audios/classical.wav');
audio=audioread('audios/hiphop.wav');
%%audio=audioread('audios/rock');
audio = sum(audio, 2) / size(audio, 2); %stereo to mono
%%%parameters

x=audio;

numharmonics=4; %number of processed harmonics
Fcc=120; % cut-off frequency (Hz)
%GdB=8; % NLD gain for transient calibration
thresh_level=-70; % magnitude threshold for peak search
inharmonicity_tolerance=0.05;% tolerance parameter for harmonic detection: fn=[n*f0*(1+tolerance) n*f0*(1-tolerance)]
alpha_low=7; % lower weighting limit in dB/oct
alpha_high=2; % upper weighting limit in dB/oct
freq_window_hz=55; % size of the frequency window for harmonic enhancement
extra_harmonic_gain=3; %Extra harmonic gain for tonal calibration

%% Splitting and resampling
%LPF
Fc=2000;
N11=2000;
SidelobeAtten=60;  % window side lobe atten
win=chebwin(N11+1,SidelobeAtten);  %chebychev window
b=fir1(N11, Fc/(Fs/2), 'low', win, 'scale');
LPF=dfilt.dffir(b);
xLL=filter(LPF,x);

%HPF
Fc=2000;
N11=2000;
SidelobeAtten=60;  % window side lobe atten
win=chebwin(N11+1,SidelobeAtten);  %chebychev window
b=fir1(N11, Fc/(Fs/2), 'high', win, 'scale');
HPF=dfilt.dffir(b);
xh=filter(HPF,x);

%Resample
Fsreduced=4096;
xL=resample(xLL,Fsreduced,Fs);

%%Fuzzy separation
s_win=512; %window size
n1=64; %hop 

S = 1;
[X,Xs,Xt,Xn,Rs,Rt,Rn] = decomposeSTN(xL,S,s_win,n1,Fsreduced);

winresynth = hann(s_win);    % Window 
yt=istft(Xt,n1,winresynth,winresynth);
yn=istft(Xn,n1,winresynth,winresynth);
ys=istft(Xs,n1,winresynth,winresynth);

%% Transient processing
%LPF and HPF  splitting
Nt=500;  %filter order
SidelobeAttent=40;  % window side lobe atten
wint=chebwin(Nt+1,SidelobeAttent);  %chebychev window
b_low=fir1(Nt, Fcc/(Fsreduced/2), 'low', wint, 'scale');
b_high=fir1(Nt, Fcc/(Fsreduced/2), 'high', wint, 'scale');
LPF=dfilt.dffir(b_low);
HPF=dfilt.dffir(b_high);
yt_low=filter(LPF,yt);
yt_high=filter(HPF,yt);

%NLD
H_hwr=[0.0278 0.5 1.8042 0 -6.0133 0 12.2658 0 -11.8891 0 4.3149];
H_fexp1=[0 1.4847 0 -1.4449 0 2.4713 0 -2.4234 0 0.9158 0];
h=H_hwr;
for n=1:11
	if n==1
		yt_low_proc=h(n);
	else
		yt_low_proc=yt_low_proc+h(n)*(yt_low.^(n-1));
	end
end
h=H_fexp1;
for n=1:11
	if n==1
		yt_low_proc_2=h(n);
	else
		yt_low_proc_2=yt_low_proc_2+h(n)*(yt_low_proc.^(n-1));
	end
end
%BPF
NBPF1     = 200;  % Order
NBPF2     = 100;  % Order
Fstop = 2*Fcc/4;   % Stopband Frequency
Fpass = Fcc+ Fcc/4;  % Passband Frequency
Fpass2=2*Fcc;
Fstop2=5*Fcc;
Wstop = 1;    % Stopband Weight
Wstop2 = 0.001;    % Stopband Weight
Wpass = 1;    % Passband Weight
Wpass2 = 1;    % Passband Weight
dens  = 20;   % Density Factor
b  = firpm(NBPF1, [0 Fstop Fpass Fsreduced/2]/(Fsreduced/2), [0 0 1 1], [Wstop Wpass], ...
           {dens});
Hd = dfilt.dffir(b);
yt_low_proc_bpf=filter(Hd,yt_low_proc_2);
b  = firpm(NBPF2, [0 Fpass2 Fstop2 Fsreduced/2]/(Fsreduced/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
Hd = dfilt.dffir(b);
yt_low_proc_bpf=filter(Hd,yt_low_proc_bpf);

%delay adjustment
delay_low=round((NBPF1+NBPF2)/2);
yt_high=[zeros(delay_low,1); yt_high];
if length(yt_high)>length(yt_low_proc_bpf)
	yt_high=yt_high(1:length(yt_low_proc_bpf));
else
	yt_low_proc_bpf=yt_low_proc_bpf(1:length(yt_high));
end

%gain calculation
N_gain=512;
R_gain=64;
pend=length(yt_low);
yt_low_padded=[yt_low; zeros(1,N_gain)'];
yt_low_proc_bpf_padded=[yt_low_proc_bpf; zeros(1,N_gain)'];
pin=0;
iii=0;
aWeight = weightingFilter('K-weighting','SampleRate',4096);
while pin<pend
	iii=iii+1;
	%chunked_yt_low(i,:)=yt_low(pin+1:pin+N_gain);
	grain_low=yt_low_padded(pin+1:pin+N_gain);
	filtered_grain_low=aWeight(grain_low);
	power_low=sum(filtered_grain_low.^2)/N_gain;
	db_low(iii)=10*log10(power_low);

	grain_proc=yt_low_proc_bpf_padded(pin+1:pin+N_gain);
	filtered_grain_proc=aWeight(grain_proc);
	power_proc=sum(filtered_grain_proc.^2)/N_gain;
	db_proc(iii)=10*log10(power_proc);

	pin=pin+R_gain;
end

difference=db_low-db_proc;

db_low_rs=resample(db_low,Fsreduced,Fsreduced/R_gain);
db_proc_rs=resample(db_proc,Fsreduced,Fsreduced/R_gain);
difference_rs=resample(difference,Fsreduced,Fsreduced/R_gain);
db_low_rs=[zeros(1,(N_gain+R_gain)/2) db_low_rs];
db_proc_rs=[zeros(1,(N_gain+R_gain)/2) db_proc_rs];
difference_rs=[zeros(1,(N_gain+R_gain)/2) difference_rs];
difference_rs(db_low_rs<thresh_level)=0;
difference_rs=smooth(difference_rs,N_gain/2);
difference_rs=difference_rs(1:length(yt_low_proc_bpf));

%Apply gain

G=10.^(difference_rs/20);
%G=G/2;
yt_low_proc_gain=G.*yt_low_proc_bpf;

%Reconstruction
yt_proc=yt_high+yt_low_proc_gain;

delay_transients=Nt/2 +delay_low;

%%Harmonic detection variables

%initialisation
[L1, nframes]=size(X);

f0found=zeros(1,nframes);
detected_f0=zeros(3,nframes);
detected_f0_values=zeros(3,nframes);
detected_harmonics=zeros(numharmonics,nframes); 
detected_harmonics_values=zeros(numharmonics,nframes);
fixed_low=zeros(L1,nframes).*(-Inf);
fixed_high=zeros(L1,nframes).*(+Inf);
accumphases=zeros(L1,1);
YL=zeros(L1,nframes);

f0max=(Fcc*s_win/Fsreduced);
f0min=f0max/4;

%%Harmonic enhancement variables

freqwindowsize=round(freq_window_hz*s_win/Fsreduced); 
if mod(freqwindowsize,2)==0
        freqwindowsize=freqwindowsize+1;
end
freqwindow=tukeywin(freqwindowsize,0.75);
freqwindow(freqwindow==0)=1e-3;



%%Harmonic weighting variables

target_weights_timbre=zeros(L1,nframes);
numBands=7;
range=[Fcc/4, Fcc*numharmonics];
[FB, cf, il, ih]=v_filtbankm(numBands,s_win,Fsreduced,range(1),range(2),'bc');
FBfull=full(FB);
FBB=zeros(floor(s_win/2 )+1, numBands);
FBB(il:ih,:)=FBfull';
FBB=FBB./sum(FBB,1);
cfbins=round(cf*s_win/Fsreduced) +1;

%%Tonal processing
for n=1:nframes
        if mod(n,100)==0
                n/nframes
        end
        f=Xs(:,n)./s_win;
        r= abs(f);	% magnitude

        bark=FBB'*r;

        [exact, exact_peak]= get_peaks(r,Fsreduced, thresh_level, Fcc/4);
        xq=cfbins(1):cfbins(end);
        interpbark=zeros(L1,1);
        interpbark(xq)=interp1(cfbins,bark,xq,'pchip');
        interpbark=20*log10(interpbark);

        locations_fundamental_all=exact(exact<f0max);

        peak_fundamental_all=exact_peak(exact<f0max);
        if ~isempty(locations_fundamental_all)

                [maxim, loc]=max(peak_fundamental_all);
                maximloc=locations_fundamental_all(loc);
                alimit=maximloc*inharmonicity_tolerance;
                borders=[1/3 1/2 1]*f0max;

                if maximloc/3 > f0min
                        f0cmin=(maximloc-alimit)/3;
                        f0cmax=(maximloc+alimit)/3;
                        [a,s]= min(abs(exact-maximloc/3));
                        if exact(s)>f0cmin & exact(s)<f0cmax 
                                a=exact(s)<=borders;
                                b=find(a,1,'first');
                                detected_f0(b,n)=exact(s);
                                detected_f0_values(b,n)=exact_peak(s);
                                f0found(n)=b;
                        end
                end

                if f0found(n)==0 & maximloc/2> f0min
                        f0cmin=(maximloc-alimit)/2;
                        f0cmax=(maximloc+alimit)/2;
                        [a,s]= min(abs(exact-maximloc/2));
                        if exact(s)>f0cmin & exact(s)<f0cmax
                                a=exact(s)<=borders;
                                b=find(a,1,'first');
                                detected_f0(b,n)=exact(s);
                                detected_f0_values(b,n)=exact_peak(s);
                                f0found(n)=b;
                        end
                end

                if f0found(n)==0
                        a=maximloc<=borders;
                        b=find(a,1,'first');
                        detected_f0(b,n)=maximloc;
                        detected_f0_values(b,n)=maxim;
                        f0found(n)=b;
                end

                for h=1:numharmonics 

                        pos=detected_f0(f0found(n),n)*(h+4-f0found(n));

                        [a,b]= min(abs(exact-pos)); alimit=pos*inharmonicity_tolerance;
                        if a<alimit
                                detected_harmonics(h,n)=exact(b);
                                detected_harmonics_values(h,n)=exact_peak(b); %value harm
                        else
                                detected_harmonics(h,n)=0;
                                detected_harmonics_values(h,n)=0; %value harm
                        end
                end

                interpbark=interpbark+(maxim+10^(extra_harmonic_gain/20)-interpbark(round(maximloc+1)));
        else
                interpbark=interpbark*(-Inf);
        end

        start_value=interpbark(floor(f0max)+1);
        fixed_low(1:floor(f0max))=-Inf;
        fixed_high(1:floor(f0max))=+Inf;
        fixed_low(floor(f0max)+1:L1,n)=start_value-alpha_low*((0:L1-floor(f0max)-1)/floor(f0max));
        fixed_high(floor(f0max)+1:L1,n)=start_value-alpha_high*((0:L1-floor(f0max)-1)/floor(f0max));

        target_weights_timbre(:,n)=interpbark;


        fsynth=f;
        shiftleft=zeros(numharmonics,1);
        shiftright=zeros(numharmonics,1);

        new_freq_bin=zeros(numharmonics,1);
        sel_weight=zeros(numharmonics,1);

        if f0found(n)~=0
                fundamental_exact=detected_f0(f0found(n),n);
                fundamental_bin=floor(fundamental_exact)+1;

                delta=floor(freqwindowsize/2);
                left=fundamental_bin-delta;
                right=fundamental_bin+delta;
                region=f(left:right);
                regionr=abs(region);
                phregion=angle(region);
                for nh=1:numharmonics
                        beta=nh+(4-f0found(n));
                        if detected_harmonics(nh,n)==0 
                                newfreq=beta*fundamental_exact;
                                binshift=newfreq-fundamental_exact;
                                shift=(binshift)*2*pi/s_win;

                                orderfracfilter=4;
                                d = fdesign.fracdelay(binshift-floor(binshift),'N',orderfracfilter);
                                secondOrderFrac = design(d,'lagrange','filterstructure','farrowfd');

                                shiftleft(nh)=left+floor(binshift);
                                shiftright(nh)=right+floor(binshift);

                                phregion=unwrap(phregion);
                                phregionpad=[phregion; zeros(orderfracfilter/2,1)];
                                phregionfilt=filter(secondOrderFrac, phregionpad);
                                phregionfilt=phregionfilt(orderfracfilter/2 +1:length(phregionfilt));
                                p0=accumphases(round(newfreq));
                                pu=p0+(n1)*shift;
                                phregionx=phregionfilt+pu;
                                accumphases(shiftleft(nh):shiftright(nh))=ones(length(region),1)*pu;

                                timbre_weight=target_weights_timbre(round(newfreq)+1,n);
                                new_freq_bin(nh)=round(newfreq)+1;
                                sel_weight(nh)=timbre_weight;

                                low_weight=fixed_low(round(newfreq)+1,n);
                                high_weight=fixed_high(round(newfreq)+1,n);
                                if timbre_weight<low_weight
                                        sel_weight(nh)=low_weight;
                                elseif timbre_weight> high_weight
                                        sel_weight(nh)=high_weight;
                                end
                                if sel_weight(nh)> thresh_level
                                        regionr2=regionr.*(10^((sel_weight(nh)-detected_f0_values(f0found(n),n))/20));
                                        regionrpad=[regionr2; zeros(orderfracfilter/2,1)];
                                        regionrfilt=filter(secondOrderFrac, regionrpad);
                                        regionrfilt=regionrfilt(orderfracfilter/2 +1:length(regionrfilt));

                                        fsynth(shiftleft(nh):shiftright(nh))=(abs(fsynth(shiftleft(nh):shiftright(nh)))+((regionrfilt-abs(fsynth(shiftleft(nh):shiftright(nh)))).*freqwindow)).*exp(i*phregionx);
                                end
                        elseif detected_harmonics(nh,n)~=0 
                                harmonic=detected_harmonics(nh,n);
                                harmonic_bin=round(harmonic)+1;
                                shiftleft(nh)=harmonic_bin-delta;
                                shiftright(nh)=harmonic_bin+delta;
                                harmonicregion_r=abs(f(shiftleft(nh):shiftright(nh)));
                                harmonicregion_phi=angle(f(shiftleft(nh):shiftright(nh)));

                                timbre_weight=target_weights_timbre(round(detected_harmonics(nh,n))+1,n);
                                new_freq_bin(nh)=round(detected_harmonics(nh,n))+1;
                                sel_weight(nh)=timbre_weight;

                                low_weight=fixed_low(round(detected_harmonics(nh,n))+1,n);
                                high_weight=fixed_high(round(detected_harmonics(nh,n))+1,n);
                                if timbre_weight<low_weight
                                        sel_weight(nh)=low_weight;
                                elseif timbre_weight> high_weight
                                        sel_weight(nh)=high_weight;
                                end
                                if detected_harmonics_values(nh,n)<sel_weight(nh) & sel_weight(nh)>thresh_level
                                        harmonicregion_r_gain=harmonicregion_r.*10^((sel_weight(nh)-detected_harmonics_values(nh,n))/20);
                                        fsynth(shiftleft(nh):shiftright(nh))=(harmonicregion_r +(harmonicregion_r_gain-harmonicregion_r).*freqwindow).*exp(i*harmonicregion_phi);
                                end
                                accumphases(shiftleft(nh):shiftright(nh))=harmonicregion_phi;
                        end
                end
                a=[];
                b=1:length(f);
                for nh=1:numharmonics
                        a=[a shiftleft(nh):shiftright(nh)];
                end
                b(a)=[];
                accumphases(b)=0;
        end

        YL(:,n)=fsynth;

end
yL=istft(YL,n1,winresynth,winresynth).*s_win;

%% transient tonal reconstrction
yL=[zeros(delay_transients,1); yL];
yn=[zeros(delay_transients,1); yn];
lengths=[length(yL) length(yn) length(yt_proc)];
[a,b]=min(lengths);
yL=yL(1:a);
yn=yn(1:a);
yt_proc=yt_proc(1:a);

yLst=yL+yn+yt_proc;

%%resample original Fs
yLst=resample(yLst,Fs,Fsreduced);


%%apply delay to HPF input
delay=ceil(s_win*Fs/Fsreduced+delay_transients*Fs/Fsreduced);
xh=[zeros(delay,1);xh];


if length(xh)>length(yLst)
        y=xh(1:length(yLst))+yLst;
else
        y=xh+yLst(1:length(xh));
end

%%Loudspeeaker simulation filter
xfilt=HPF_final(x,Fcc,Fs);
yfilt=HPF_final(y,Fcc,Fs);
y_filt_low=LPF_final(x,Fcc,Fs); 

delay_end=delay+N11/2+3000/2;
yfilt=yfilt(delay_end:end);%resulting signal filtered
delay_end_2=3000/2;
xfilt=xfilt(delay_end_2:length(yfilt)+delay_end_2-1); %high-passed version
y_filt_low=y_filt_low(delay_end_2:length(yfilt)+delay_end_2-1);
y_darre=y_filt_low+yfilt; % resulting signal with low components

if length(x) > length(y_darre)
    x2=x(1:length(y_darre)); %original signal
end

