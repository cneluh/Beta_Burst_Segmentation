function xf = fold(x,xl,xr,bu,bl)
% xf = fold(x,xl,xr,,bu,bl)
% x,xl,xr,bu,br Signal, left, right overlapping regions
% bu, bl Upper and lover part of smooth rising function
% xf  the folded data

	m = length(bl);
	n = length(x);
    
    nl=length(xl);
         
    front = 1:m;
	back  = n:-1:(n-m+1);
    xl_back  = nl:-1:(nl-m+1);
    
    x(front)  = x(front).*bu+ xl(xl_back).*bl;
    x(back)  = x(back).*bu-xr(front).*bl;
    xf=x;

% ALCB Version 1.0
% Created, 29.10.2003
% e-mail firat@eemb.cu.edu.tr
% N.Firat INCE

    
