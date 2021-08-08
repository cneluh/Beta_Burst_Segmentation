function extra = edge_fold(which,xc,bp,bm)
% edgefold -- Perform folding projection with (+,-) polarity at EDGES
%  Usage
%    extra = edgefold(which,xc,bp,bm)
%  Inputs
%    which  string, 'left'/'right', indicating which edge we are at
%    xc     unfolded data, center window
%    bp     interior of window
%    bm     exterior of window
%  Outputs
%    extra  pseudo-left/right packet
%
%  Description
%    The result should be used as either left or right packet in
%    fold to ensure exact reconstruction at edges.
%
%  See Also
%    fold, unfold, CPSynthesis, CPImpulse
%

	n = length(xc);
	m = length(bp);
	back  = n:-1:(n-m+1);
	front = 1:m;
	extra = xc.*0; 
%
	if strcmp(which,'left'),
		extra(back) = xc(front) .* (1-bp)./(bm);  % Calculates the edge data on the first and last segments
    else                                                      % In the fold function this data helps the start and end points remains the same as input  
		extra(front) = -xc(back) .* (1-bp)./bm;
	end	

%
%  Copyright (c) 1994. David L. Donoho
%
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
