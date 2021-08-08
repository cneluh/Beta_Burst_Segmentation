function [smb]=smbell(left_edge,center,right_edge)
%[smb]=smbell(left_edge,center,right_edge)
%trns:length of transisiton =2.^n
%l: length of segment 1
%Rising Function r(t)= sin[pi/4*(1+sin(pi/2*t)];
%The interval for left part is[-1 1] and for right part [1 3]

if left_edge==right_edge
    smb=smooth_bell(left_edge,center);
    
elseif left_edge==0
    le=right_edge/2;
    t=[-le:le-1]./le;
    br=sin(pi/4*(1+sin(pi/2*(t+2))));
    smb=[ones(1,center) br];
    
elseif right_edge==0
    le=left_edge/2;
    t=[-le:le-1]./le;
    bl=sin(pi/4*(1+sin(pi/2*t)));
    smb=[bl ones(1,center)];
end