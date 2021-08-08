function [ts,tsm,segind,et_seg,et_total]=lcpSeg_nondyadic(x,windw,fs,typE)
%[ts,tsm]=lcpSeg_nondyadic(x,d,bu,bl,fs)
% Returns the adaptive time segments (ts and tsm) estimated by Local Cosine Packets
% via entropy minimization/best basis algorithm
% % Use coefficient averaging if x is a matrix (multiple trials)
% x : Signal to be segmented
% windw : Smallest window size of a segment
% lcp coefficiets are used to calcute entropy and chose segmentation
% typE = 1 : FDCTIV
%		 0 : FFT

[M,N] = size(x);
segLength =windw*fs; % window size in sample size
[bl,bu]=smooth_edge(segLength/2); % Smooth Orthagonal Edges
fprintf('Smooth Edge Size=%d\n',segLength/2);
ts=[];
tsm=[];
nOfSeg = N/segLength; %this gives me the the total number of smaller segments
et = zeros(size(x,2),3);

% check legth of data and window size
  if mod(nOfSeg,1) ~= 0
      error('Length of data cannot be divided evenly with window size');
  end
  
% LCP of entire signal 
lcp_dcx=0;
for xx = 1:1:size(x,1) 
  x_=x(xx,:);
  xl=edge_fold('left',x_,bu,bl);
  xr=edge_fold('right',x_,bu,bl);
  x1=fold(x_,xl,xr,bu,bl);
  lcp_dcx= lcp_dcx + fdctIV(x1).^2; % DCT of entire signal
end
  lcp_dcx1= (lcp_dcx)/sum(lcp_dcx); % normalizing
  ss = norm(sqrt(lcp_dcx1));
  et_total = entropy(sqrt(lcp_dcx1)/ss);

% Start calculation
et_seg=[];
segind = zeros(nOfSeg-1,4);% saving the index of each fixed window
for i=0:1:nOfSeg-2 % fixed window
        
   
dcl=0;
dcm=0;
dcf=0;
for xx = 1:1:size(x,1)
    xs = x(xx,:);% choosing which shifted data to compute first
    
    if i==0 % gives intial starting points for first child 
    start_ = i*segLength+1;
    middle_child=(i+1)*segLength;
    end
    
    middle_fixed =(i+1)*segLength+1; % would starting point of the fixed window
    end_= ((i+2)*segLength); % would be the end point of fixed window
    
 
 % extracting the segments
 xleft = xs(start_:middle_child);
 xfixed = xs(middle_fixed:end_);
 xmother = xs(start_:end_);
 mother_length = length(xmother);  
 
 %Left Edge
    if start_ ==1
         xl_left=edge_fold('left',xleft,bu,bl);
    else
         xl_left=xs((i-1)*segLength+1:i*segLength);
    end
    
    %Since length of mother varies, we must change the size of the left edges
    %accordingly
    
    if start_ ==1 || mother_length > start_
        xl_mother = edge_fold('left',xmother,bu,bl);
    else
        xl_mother=xs((start_-length(xmother)):(start_-1));
    end
    
    xl_fixed=xs(i*segLength+1:(i+1)*segLength);% left edge of fixed

 %Right Edge
    if i+1<nOfSeg-1
         xr_fixed=xs((i+2)*segLength+1:(i+3)*segLength);
    else
         xr_fixed=edge_fold('right',xfixed,bu,bl);
    end
    
    %Since length of mother varies, we must change the size of the right edges
    %accordingly
    
    if (end_ + mother_length)<= length(xs)
        xr_mother= xs(end_+1:(end_+mother_length));
    else
        xr_mother=edge_fold('right',xmother,bu,bl);
    end

  xr_left=xs((i+1)*segLength+1:(i+2)*segLength); 
   
 %Construct the segments w/ appropiate edges
    xfixed=fold(xfixed,xl_fixed,xr_fixed,bu,bl);
    xleft=fold(xleft,xl_left,xr_left,bu,bl);
    xmother=fold(xmother,xl_mother,xr_mother,bu,bl);

    if typE==1
       
    dcfixed=fdctIV(xfixed).^2;
    dcleft=fdctIV(xleft).^2;
    dcmother=fdctIV(xmother).^2;

    dcl = dcl+dcleft;
    dcf = dcf+dcfixed;
    dcm= dcm+dcmother;
     
    end

end

if typE ~=1
    %computing fft using 2*seglength
    dcfixed=abs(fft(xfixed,2*length(xfixed))).^2;   % Fourier Extended to Dyadic Length
    dcfixed = dcfixed(1:length(xfixed))./length(xfixed);

    dcleft=abs(fft(xleft,2*length(xleft))).^2; 
    dcleft = dcleft(1:length(xleft))./length(xleft);

    dcmother=abs(fft(xmother,2*length(xmother))).^2; 
    dcmother=dcmother(1:mother_length)./mother_length;
    
    dcl = dcleft;
    dcf = dcfixed;
    dcm= dcmother;
end
    

%second method
ss= sum(dcm); % sum of the mother
lcp_fixed= dcf/ss;
lcp_left= dcl/ss;
lcp_mother= dcm/ss;

 if typE ==1
    [et]= entropy_nondyadic_tree(sqrt(lcp_left),sqrt(lcp_fixed),sqrt(lcp_mother));
 else
    et = zeros(1,3); % preallocating entropies
    et(1) = entropy(sqrt(lcp_left));
    et(2) = entropy(sqrt(lcp_fixed));
    et(3)= entropy(sqrt(lcp_mother));
 end


% Best basis algorithm between mother and children

if et(1)+et(2) < et(3)
       % keep child
       ts = [ts; (start_+length(xleft)-1)/fs]; % would the length of the segmented data in seconds
       tsm = [tsm ; start_ , middle_child]; % would be the index where the data is segmented
       %et_seg= [et_seg,et(1)+et(2)];
       et_seg= [et_seg,et(1)];
       
       start_= middle_fixed;
       middle_child = end_;
            
   else % if entropy of the mother is less than the children combined we keep the mother
     % start would be same
     middle_child = end_;
          
   end
   
   % save index for last segment of fixed and mother
   if et(1)+et(2) < et(3) && i == nOfSeg-2
       ts = [ts ; (start_+length(xfixed)-1)/fs];
       tsm = [tsm ; middle_fixed, end_];
        %et_seg= [et_seg, et(1)+et(2)];
        et_seg= [et_seg, et(2)];
   elseif et(1)+et(2) > et(3) && i == nOfSeg-2
        ts = [ts ; (start_+mother_length-1)/fs];
       tsm = [tsm ; start_, end_];
       et_seg= [et_seg, et(3)];
       
   end

end

ts=ts';


