function [bl, bu]=smooth_edge(n)
% function [bl bu]=smooth_edge(n)
% Smooth edges for folding
% n=length of smooth rising function

t=[0:n-1]./n; % [0 1] interval
bu=sin(pi/4*(1+sin(pi/2*t))); % Calculates rising function left upper part
bl=sqrt(1-bu.^2); % Orhagonal lover part.

% ALCB Version 1.0
% Created, 29.10.2003
% e-mail firat.ince@gmail.com
%  N.Firat INCE
    
