function [tf,t,f,lsp]=CreateTFMap2(cn,ts,fs,tyP,SmoothSize,edge_,NFFT)
%[tf,t,f]=CreateTFMap2(cn,ts,fs,tyP,SmoothSize,edge_)
%cn= Input signal
%ts= time segments (can be computed with lcpSeg_nondyadic)
%fs= sampling Frequency
%Inputs especifically for LocalSmoothPackets2:
%tyP=   1: DCT-IV
%       2: FFT
%edge_= 1:Folded data
%       2:Masked data

tsF=ts;
tsL=ts;% to inputs original ts into LocalSmoothPackets2

tsF=(tsF*fs);
[m,~]=size(cn);
lts=length(tsF);
MaxTimeSeg= max([tsF(1) diff(tsF)]);  % Find number of FFT points = Maximun Time Segment
MinTimeSeg= min([tsF(1) diff(tsF) 256]);  % So that all TF maps have sime size (Time)


st=0;
f=(0:NFFT/2-1)./(NFFT/2)*fs/2;% frequency
tsF=round([0 tsF]); % Reshape the time segment with zero time index

tf=zeros(NFFT/2,m/MinTimeSeg); % Generates the TF map having 'lts' colums (creates matrix of zeros for Tf map, where the values of each frequency(row) and minute(colunm)
t=[0:(m/MinTimeSeg-1)]*MinTimeSeg/fs; % time vector

[lsp,~]=LocalSmoothPackets2(cn,tsL,fs,tyP,SmoothSize,edge_,NFFT); 

for i=2:lts+1
    
    SegLen=tsF(i)-tsF(i-1);

    coef=lsp(:,i-1);
    coef=coef(1:NFFT/2,:)/(SegLen);
    coef=mean(coef,2); % for generated data

    %compute coefficient to remap it to fit matrix
    tcoef=coef;
    tcoef=repmat(tcoef,1,SegLen/MinTimeSeg);
    
    tf(:,[st/MinTimeSeg+1:tsF(i)/MinTimeSeg])=tcoef;
    st=tsF(i);
end

