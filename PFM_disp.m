%% 06-05-2021
%% WHEN PHI PART DISABLED WORKS AS ELASTIC PROBLEM
clear all
close all
clc

%declare global variables here
global E nu P C
global node element numnode numelem ndof tdof
global Gc Lc ks

%------------------------------- Meshing----------------------------------
node=load('node.txt');
element=load('element.txt');
crack = [0.0 0.0; 0.5 0.0] ;
lsm=getLevelSetFunc(crack);

numnode = size(node,1) ;
numelem = size(element,1) ;
elemnode=ones(numelem,1)*4;

%% Material properties
E = 210000; nu = 0.3; P = 1; R=0;
ks=1e-11; % stability term
Gc=2.7;Lc=.7e-1;
C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];

% loop over number of steps
% ndof = 2 ;
% ndofPF=3; %x1,y1,x2,y2,...xn,yn,phi1,phi2...,phin
ndof = 3 ;
tdof = numnode*ndof;
qvec=zeros(tdof,1);
fvec = zeros(tdof,1) ; % External. Zero in disp control

botNodes = reshape(find(node(:,2)==-0.5),[],1);
topNodes = reshape(find(node(:,2)==0.5),[],1);

bcdof = []; bcval = [];
bcreactX =[]
bcreactY =[]

% Y and X on bottom
for i = 1:length(botNodes)
    bcdof = [bcdof 2*botNodes(i)-1 2*botNodes(i)];
    bcval = [bcval 0 0];
    bcreactX = [bcreactX 2*botNodes(i)-1];
    bcreactY = [bcreactY 2*botNodes(i)];
end
%bcreact=bcreactY;
bcreact=[bcreactX bcreactY];
rvec = zeros(tdof,1) ;
%% Begin load steps
ipas = 1 ;
Lambda = 10; % Displacement increment
Papp = 1;
ipas = 1
xdof=(1:2:2*numnode);
ydof=(2:2:2*numnode);
udof=[1:2*numnode];
tu=length(udof);
phidof=(tu+1:tdof);
tphi=length(phidof);
actDOF=setdiff([1:tdof],bcreact); % active U and PHI. Neeed for residual
actDOF_U=setdiff(udof,bcreact); % active U

tuact=length(actDOF_U);
tphi=length(phidof);

qvec=zeros(tdof,1); % 3*numnode
qu=zeros(tu,1); % 2*numnode
qphi=zeros(tphi,1);

fvec=zeros(tdof,1);% 0 in displacement control
fu = zeros(tu,1); 
fphi=zeros(tphi,1);

ress=[];

%%%%%%%%%%%% DISPLACEMENT COMTROL
dispdof=2*topNodes;
PFMPHI=[];
Pmax=0;
Umax=0;
disp_y=0;
phi_y=0;
force_y=0;
Uapp=0.0;
LambdaD=0.001;
while Umax<0.015
    Uapp=Uapp+LambdaD;
    disp([num2str(toc),'  Force step  ',num2str(ipas)]) ;
    iter=1;
    err=1;
    errU=1;
    errPHI=1;
    QPHI=zeros(tphi,1);
    QU=zeros(tuact,1);
    while err>1e-10
        qvec1=qvec; % qvec input to this iteration step
        qvec2=qvec; % qvec output from this iteration step.       
        KTOT=zeros(tdof,tdof);
        RTOT=zeros(tdof,1);
        
        %%%%%%SOLVE FOR PHI===================================
        [KTOT,RTOT]=getPFMRes(qvec1,KTOT,RTOT);%        
        QPHI=QPHI-KTOT(phidof,phidof)\RTOT(phidof);
        qvec2(phidof,1)=QPHI; % PHI solution updated in qvec2
        qvec=qvec2; % Updated
        
        %%%%%%SOLVE FOR U===================================
        [KTOT,RTOT]=getElasticRes(qvec2,KTOT,RTOT);        
        for i=1:length(dispdof)
            c=dispdof(i); %current dof
            loadvalue(i)=Uapp;
        end
        
        bcdofmix=[bcdof reshape(dispdof,1,[])]; % FixedDOF+DiscontrolDOF
        bcvalmix=[bcval loadvalue];
        [kmod,fmod]=feaplyc2(KTOT,RTOT,bcdofmix,bcvalmix);
        QU=kmod(actDOF_U,actDOF_U)\fmod(actDOF_U);    
        qvec2(actDOF_U,1)=QU;

        res=qvec2-qvec1;
        err= norm(res);    
        disp(['Step No: ' num2str(ipas) 'Iteration No: ' num2str(iter) ' Error:Phi '  num2str(err)]);
        iter=iter+1;
        qvec=qvec2;
        itermax=10;
        if iter>itermax
            break
        end
    end
    ress=[ress;[ipas iter errU errPHI]];
    fbulk=KTOT*qvec;
    stdux = qvec(1:2:2*numnode);
    stduy = qvec(2:2:2*numnode) ;
    stdphi = qvec(2*numnode+1:tdof) ;
    disp_y(ipas) = Uapp;
    phi_y(ipas) = max(stdphi);
    fbulk=KTOT*qvec;
    force_y(ipas)= sum(fbulk(dispdof));
    Pmax(ipas)=force_y(end);
    Umax(ipas)=disp_y(end);
    tt(ipas)=trace(KTOT);
    PFMPHI=[PFMPHI stdphi];
    ipas=ipas+1;
end % End Force increment

figure(3)
%% plot P-delta
load('Emilio.txt');
plot([0 disp_y(1:end)],[0 force_y(1:end)],'r-o');hold on
plot(Emilio(:,1),Emilio(:,2),'b-o');
plot([0 0.005],[0 719.46],'k--') % Elastic solution
legend('Present', 'Elmilio', 'Elastic')