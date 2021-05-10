function [node,element]=separationCrackNode_lsm_func_v1(L,W,xCr,node,element,nxdiv,lsm,clen) 
anode=clen/(W/(nxdiv-1)) % additional nodes
numnode=size(node,1);
k=1;
node_ex=[];  % extra nodes
node_num=[]; % node index
for i=1:numnode
    if lsm(i,1)==0 && lsm(i,2)<0 % nodes to left of left crack. Needs to be duplicated
%         node(numnode+k,:)=node(i,:)
        node=[node;node(i,:)];
        node_ex=[node_ex;node(i,:)];
        node_num=[node_num;i]; % number of original nodes
        k=k+1;
    end
end
node_dup=[1:size(node_num,1)]'+numnode; % numbering of duplicate nodes

numelem=size(element,1);
k=1;
for i=1:numelem
    nix=element(i,:);
    lsm_cg=sum(lsm(nix,:))/3;
    comm=intersect(nix,node_num); % find in which element the nodes are attached
    cl=length(comm);
    if cl~=0 && lsm_cg(1)>0 && lsm_cg(2)<0 % has crack nodes,left&above crack
        for ix=1:cl
            a=find(node_num==comm(ix)); % index of common node that is duplicated
            b=find(element(i,:)==comm(ix));
            element(i,b)=node_dup(a);
        end
    end
end