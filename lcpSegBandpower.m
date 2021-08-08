function [SegPower,f] = lcpSegBandpower(xx,ts,fs,tYp,Smoothsize,edge_,resltn,freqRange)
% [lspm,tsm]=lcpbandpower(xx,ts,fs,tyP,SmoothSize,edge_,resltn)
% ts= is the adaptive segmentation returned from  LDB or BestBase
% fs= Sampling Frequency
% tyP =    1 : DCT-IV
%          2 : Fourier Smooth
% SmoothSize= the length of the overlapping part of smooth bell/Folded data
% lspm= Local Packets of adaptive time segmentation
% tsm= Time segmentation matrix Nx2
% edge_ =   1: Folding the data
%           2: Masking the data
%resltn = frequency resolution
%freq = frequency range , exmpl. [low high] = [0 50]

NFFT = resltn*fs; % nfft points for desired freqeuncy resolution
f=(0:NFFT/2-1)./(NFFT/2)*fs/2;% frequency

[lcp,tsm]=LocalSmoothPackets2(xx,ts,fs,tYp,Smoothsize,edge_,NFFT); %Local Cosine Packets

SegPower = zeros(1,length(tsm)); % Power Of each segment
for ii = 1:length(tsm)
    
    LCP = lcp(1:NFFT/2,ii);
    if length(freqRange) ==2
        freqIndx = find( f>=freqRange(1) & f<=freqRange(2));
        SegPower(ii) = sum(LCP(freqIndx));
    elseif length(freqRange) ==1   
        freqIndx = find(f==freqRange);
        SegPower(ii) = sum(LCP(freqIndx));
    end
end
       
    
    