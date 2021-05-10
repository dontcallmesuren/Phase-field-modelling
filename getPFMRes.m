function[kbulk,rvec]=getPFMRes(qvec,KTOT,RTOT)

global node element
% Global stiffness
kbulk = KTOT;
rvec = RTOT;
global numelem elemnode node element numnode ndof tdof
for iel = 1:numelem
    clear sctr sctrB sctrU;
    sctrB=[];
    sctrU=[];
    sctrPHI=[];
    sctr = element(iel,1:4);
    nn = length(sctr) ;
    nodes = node(sctr,:);
    for k = 1 : nn
        sctrB(2*k-1) = 2*sctr(k)-1 ;
        sctrB(2*k) = 2*sctr(k) ;
        sctrB(2*nn+k)   = 2*numnode+sctr(k);
        sctrU=[sctrU 2*sctr(k)-1 2*sctr(k)];
        sctrPHI=[sctrPHI 2*numnode+sctr(k)];
    end
% For regular FEM
    IntOrder = 2 ;
    q=qvec(sctrB); % u and phi
    u=qvec(sctrU);
    phi=qvec(sctrPHI);
    [elmat,rmat] = getkbulk_res_PHI(iel,IntOrder,u,phi);
    kbulk(sctrB,sctrB) = kbulk(sctrB,sctrB) + elmat ;
    rvec(sctrB,1) = rvec(sctrB,1) + rmat ;
end