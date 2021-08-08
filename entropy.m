function e=entropy(x);
% calculates Shannon Entropy
% e=entropy(x);
x=x.^2;
e=-sum(x.*log(x+eps));