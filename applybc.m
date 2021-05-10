function [bcdof,bcval]=applybc(ndof,bcdof,bcval,topNodes,botNodes,Uapp)

for i = 1:length(topNodes)
    bcdof=[bcdof ndof*topNodes(i)];
    bcval = [bcval Uapp];
end

%for i = 1:length(botNodes)
%    bcdof=[bcdof ndof*botNodes(i)];
%    bcval = [bcval -Uapp];
%end
