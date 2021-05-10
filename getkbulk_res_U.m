function[elmat,rmat]=getkbulk_res_U(iel,IntOrder,u,phi)
global node element elemType
global BG kmat E1 E2 C C1 C2 Rinc
global stressState nu
global sigdata
global Gc Lc ks
global numnode ndof tdof
% Element stiffness

sctr = element(iel,1:4);
cg = mean(node(sctr,:));
nn = length(sctr) ;
ndof=3;

intType = 'GAUSS' ;
[W,Q] = quadrature(2) ;
elmat=zeros(3*nn,3*nn);
kuu=zeros(2*nn,2*nn);
kpp=zeros(1*nn,1*nn);

rmat=zeros(3*nn,1);
ru=zeros(2*nn,1);
rphi=zeros(nn,1);

for kk = 1:size(W,1)
    Gpt = Q(kk,:) ;
    [N,dNdxi] = lagrange_basis(Gpt) ; % element shape functions
    Gxy = N'*node(sctr,:) ;
    J0 = node(sctr,:)'*dNdxi ; % element Jacobian matrix
    invJ0 = J0\eye(2);
    dNdx  = dNdxi*invJ0;       % derivatives of N w.r.t XY
    
    %Standard B matrix
    Bu = zeros(3,2*nn) ;
    Bu(1,1:2:2*nn) = dNdx(:,1)' ;
    Bu(2,2:2:2*nn) = dNdx(:,2)' ;
    Bu(3,1:2:2*nn) = dNdx(:,2)' ;
    Bu(3,2:2:2*nn) = dNdx(:,1)' ;
    
    BPhi = zeros(2,1*nn) ;
    BPhi(1,1:nn) = dNdx(:,1)' ;
    BPhi(2,1:nn) = dNdx(:,2)' ;
    
    % Non linear terms
    phix=N'*phi; % Phi at Gauss point (x)
    delPHI=BPhi*phi; % vector [dphi/dx; dphi/dy] at Gauss point
    e0=Bu*u;
    sig0=C*e0;
    SED0=0.5*e0'*C*e0;

%    if phix>0.9999
%      gphi=(1-phix)^2+ks;
%      dgphi=-2*(1-phix);
%    else
%      gphi=(1-phix)^2;
%      dgphi=-2*(1-phix);
%    end
%    if phix<=0 % CHECK LATER
%      gphi=1;
%      dgphi=-2;
%    end

    gphi=(1-phix)^2;
    dgphi=-2*(1-phix);
    if phix>0.9
      phix=1;
      gphi=ks;
      dgphi=0;
    end    
    if phix<0
      phix=0;
      gphi=1;
      dgphi=0;
    end

    kuu=kuu+gphi*Bu'*C*Bu*W(kk)*det(J0);    
    ru=ru+gphi*Bu'*sig0*W(kk)*det(J0);
    
end
sctrPhiloc=[9:12];
sctrUloc=setdiff([1:12],sctrPhiloc);
elmat(sctrUloc,sctrUloc)=kuu;
rmat(sctrUloc,1)=ru;