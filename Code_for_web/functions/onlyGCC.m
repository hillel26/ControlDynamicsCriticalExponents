function [A] = onlyGCC(A)
%ONLYGCC Summary of this function goes here
%   Detailed explanation goes here
cb = conncomp(graph(A));
[~,indGcc] = max(histcounts(cb,1:max(cb)+1));
outgcc = (cb~=indGcc);
A(outgcc,:)=[];
A(:,outgcc)=[];
end

