function [et]= entropy_nondyadic_tree(lcp_left,lcp_fixed,lcp_mother)
% computes the Shannon Entropy for nondyadic tree

et = zeros(1,3); % preallocating entropies

ss = norm(lcp_mother); % normalizing with entire signal

et(1) = entropy(lcp_left/ss);
et(2) = entropy(lcp_fixed/ss);
et(3)= entropy(lcp_mother/ss);