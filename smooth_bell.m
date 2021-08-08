function bell=smooth_bell(trns,l)
%function smooth_bell(trns,l)
%trns:length of transisiton =2.^n
%l: length of segment 1
%Rising Function r(t)= sin[pi/4*(1+sin(pi/2*t)];
%The interval for left part is[-1 1] and for right part [1 3]
trns=trns/2;
t=[-trns:trns-1]./trns;
bl=sin(pi/4*(1+sin(pi/2*t)));
br=sin(pi/4*(1+sin(pi/2*(t+2))));
bell=[bl ones(1,l) br];
