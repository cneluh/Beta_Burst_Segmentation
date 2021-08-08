function [lspm,tsm]=LocalSmoothPackets2(cn,ts,fs,tyP,SmoothSize,edge_,NFFT)
% [lspm,tsm]=LocalSmoothPackets2(cn,ts,fs,tyP,SmoothSize,edge_)
% cn is the input signal each column is a realization and
% each row are observations
% ts= is the adaptive segmentation returned from BestBase (lcpSeg_nondyadic) or LDB 
% fs= Sampling Frequency
% tyP =    1 : DCT-IV
%          2 : Fourier Smooth
% SmoothSize= the length of the overlapping part of smooth bell/Folded data
% lspm= Local Packets of adaptive time segmentation
% tsm= Time segmentation matrix Nx2
% edge_ =   1: Folding the data
%           2: Masking the data
%NFFT= points for FFT, if not using place '~';
 
%parameters for both----------
ts=round([0 ts]*fs);
lts=length(ts);
[~,M]=size(cn);
%NFFT =fs; % FFT computed over 8 seconds

if tyP ==1
lspm1=zeros(NFFT,lts);
else
lspm1=zeros(NFFT/2,lts);
end

%Parameters for folded data padded with zeros --------
MaxSeg=max(diff(ts));
MinSeg=min(diff(ts));
tsm=zeros(lts-1,2);
[bl,bu]=smooth_edge(SmoothSize);

%folded data
for m=1:M

    for i=2:lts
        
        % Both
        SegLen=ts(i)-ts(i-1);
        tsm(i-1,:)=[ts(i-1)+1,ts(i)];
        data=cn(:,m)';
        
        % Folded data ------------------------
        if edge_ == 1
        
        xc=data(ts(i-1)+1:ts(i));

        if tyP                   % Smooth DCT-IV Base
            if ts(i-1)==0
                xl=edge_fold('left',xc,bu,bl);
            else
                xl=data(ts(i-1)-MinSeg+1:ts(i-1));
            end

            if i<lts
                xr=data(ts(i)+1:ts(i)+MinSeg);
            else
                xr=edge_fold('right',xc,bu,bl);
            end

            xc=fold(xc,xl,xr,bu,bl);        % Signal is folded before DCT-IV
        end
        
       
        elseif edge_ == 2

        % Masked Data ----------------------
        if i ==2 
            if i == lts
                xc=data(ts(i-1)+1:ts(i)); % extracts the data of a
                sbell = smbell(0,SegLen,0); % computes bell with appropiate edges
            else
                xc=data(ts(i-1)+1:ts(i)+SmoothSize);
                sbell = smbell(0,SegLen-SmoothSize,SmoothSize*2);
            end

        elseif i == lts
            xc=data(ts(i-1)+1-SmoothSize:ts(i));
            sbell = smbell(SmoothSize*2,SegLen-SmoothSize,0);
        else
            xc=data(ts(i-1)+1-SmoothSize:ts(i)+SmoothSize);
            sbell = smbell(SmoothSize*2,SegLen-(SmoothSize*2),SmoothSize*2); 
        end

        xc =xc.*(sbell);% data masked with bell
        
        end
    
%calcualtion of DCT-IV and FFT -----------  

        if tyP==1
            
            r = [xc,zeros(1,(NFFT)-length(xc))]; % WITH zero padding
            dcx=fdctIV(r).^2;          % DCT-IV
            lspm1(:,i)=dcx';
            
        elseif tyP==2
            dcx=abs(fft(xc,NFFT)).^2;   %FFT
            dcx = dcx(1:NFFT/2)./(NFFT/2); % removing the symmetric side of FFT
            lspm1(:,i)=dcx';
        end
            lspm=lspm1(:,2:end);
      
    end

end

