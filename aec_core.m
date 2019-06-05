% WebRtc AEC Core Processing By Low Complex Calculating
% Luks,Jun.5th,2019
% TODO: The Adaptive Filter Of WebRtc_AEC Is UseLess.
%       Replacing New Adaptive Filter Or Fixing Old Adaptive Filter With DT
%
clc
clear
filename = '原始信号.wav';
outfilename = ['mout_core_' filename];

[sig,fs] = audioread(filename);
rrin = sig(:,2);%far-end spk
ssin = sig(:,1);%near-end mic

mult=fs/8000;
if fs == 8000
    cohRange = 2:3;
elseif fs==16000
    cohRange = 2;
end

% Flags
AECon=1;

N = 128;  % Partition length
FFT_Factor = N;%256>>1
if fs == 8000
    mufb = 0.6;
else
    mufb = 0.5;
end
%%
if fs == 8000
    echoBandRange = ceil(300*2/fs*N):floor(1800*2/fs*N);%300Hz~1.8KHz
    
else
    echoBandRange = ceil(300*2/fs*N):floor(1800*2/fs*N);%300Hz~1.8KHz
    
end
suppState = 1;
len=length(ssin);

ercn=zeros(len,1);
zm=zeros(N,1);

NN=len;
Nb=floor(NN/N);
start=1;
xo=zeros(N,1);
do=xo;
eo=xo;

dIdxV=zeros(Nb+1,1);
SLxV=zeros(Nb+1,1);
hnlSortQV=zeros(Nb+1,1);
hnlPrefAvgV=zeros(Nb+1,1);
mutInfAvg=zeros(Nb+1,1);
hnled = zeros(N+1, 1);
weight=zeros(N+1,1);
hnlMax = zeros(N+1, 1);
hnl = zeros(N+1, 1);
overdrive = ones(1, N+1);
WFbD=ones(N+1,1);
mbuf=zeros(2*N,1);

hnlLocalMin = 1;
cohxdLocalMin = 1;
ovrd = 2;
ovrdSm = 2;
hnlMin = 1;
hnlMinCtr = 0;
hnlNewMin = 0;


wins=[0;sqrt(hanning(2*N-1))];
cohxd = zeros(N+1,1);
Sd = zeros(N+1,1);
Sx = zeros(N+1,1);
Sxd = zeros(N+1,1);
hid = waitbar(0,'Processing...');
tic
for kk=1:Nb
    pos = N * (kk-1) + start;
    
    xk = rrin(pos:pos+N-1);
    dk = ssin(pos:pos+N-1);
    
    xx = [xo;xk];
    xo = xk;
    
    dd = [do;dk];
    do = dk;
    window = wins;
    if fs == 8000
        gamma = 0.9;
    else
        gamma = 0.93;
    end
    
    tmp = fft(xx.*window)/FFT_Factor;
    xf = tmp(1:N+1);
    tmp = fft(dd.*window)/FFT_Factor;
    df = tmp(1:N+1);
    
    Sd = gamma*Sd + (1 - gamma)*real(df.*conj(df));
    Sx = gamma*Sx + (1 - gamma)*real(xf.*conj(xf));
    % coherence
    Sxd = gamma*Sxd + (1 - gamma)*xf.*conj(df);
    
    cohxd = real(Sxd.*conj(Sxd))./(Sx.*Sd + 1e-26);
    hnled = 1 - cohxd;%% 越小 回声越大
    
    hnlSortQ = mean(1 - cohxd(echoBandRange));
    
    [hnlSort2, hnlSortIdx2] = sort(hnled(echoBandRange));
    hnlQuant = 0.75;%升序 取3/4处
    hnlQuantLow = 0.5;%升序 取2/4处
    qIdx = floor(hnlQuant*length(hnlSort2));
    qIdxLow = floor(hnlQuantLow*length(hnlSort2));
    hnlPrefAvg = hnlSort2(qIdx);
    hnlPrefAvgLow = hnlSort2(qIdxLow);
    
    if hnlSortQ > 0.9
        suppState = 0;
    elseif hnlSortQ < 0.8
        suppState = 1;%其他情况 保持上一帧状态
    end
    
    if hnlSortQ < cohxdLocalMin & hnlSortQ < 0.75
        cohxdLocalMin = hnlSortQ;
    end
    
    if cohxdLocalMin == 1
        ovrd = 3;
        hnled = 1-cohxd;
        hnlPrefAvg = hnlSortQ;
        hnlPrefAvgLow = hnlSortQ;
    end
    
    if suppState == 0
        hnled = 1-cohxd;
        hnlPrefAvg = hnlSortQ;
        hnlPrefAvgLow = hnlSortQ;
    end
    
    if hnlPrefAvgLow < hnlLocalMin & hnlPrefAvgLow < 0.6
        hnlLocalMin = hnlPrefAvgLow;
        hnlMin = hnlPrefAvgLow;
        hnlNewMin = 1;
        hnlMinCtr = 0;
    end
    
    if hnlNewMin == 1
        hnlMinCtr = hnlMinCtr + 1;
    end
    if hnlMinCtr == 2
        hnlNewMin = 0;
        hnlMinCtr = 0;
        ovrd = max(log(0.01)/(log(hnlMin + 1e-26) + 1e-26), 5);
    end
    hnlLocalMin = min(hnlLocalMin + 0.0008/mult, 1);
    cohxdLocalMin = min(cohxdLocalMin + 0.0004/mult, 1);
    
    if ovrd < ovrdSm
        ovrdSm = 0.99*ovrdSm + 0.01*ovrd;%release
    else
        ovrdSm = 0.9*ovrdSm + 0.1*ovrd;%attack
    end
    
    dkEn = sum(Sd);
    
    %% Overdrive
    aggrFact = 0.3;
    wCurve = [0; aggrFact*sqrt(linspace(0,1,N))' + 0.1];
    weight = wCurve;
    
    hnled = weight.*min(hnlPrefAvg, hnled) + (1 - weight).*hnled;
    od = ovrdSm*(sqrt(linspace(0,1,N+1))' + 1);
    sshift = ones(N+1,1);
    hnled = hnled.^(od.*sshift);
    hnl = hnled;
    
    if(AECon==1)
        df = df.*(hnl);
    end
    
    Fmix = df;
    % Overlap and add in time domain for smoothness
    tmp = [Fmix ; flipud(conj(Fmix(2:N)))];
    mixw = wins.*real(ifft(tmp))*FFT_Factor;
    mola  = mbuf(end-N+1:end) + mixw(1:N);
    mbuf = mixw;
    ercn(pos:pos+N-1) = mola;
    
    if mod(kk, floor(Nb/100)) == 0
        str=['Processing(' num2str(round(toc)) 's):' num2str(floor(kk/Nb*100)) '%'];
        waitbar(kk/Nb,hid,str);
    end
end
close(hid);
audiowrite(outfilename,ercn,fs);