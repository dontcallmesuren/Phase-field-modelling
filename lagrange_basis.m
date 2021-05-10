function [Nv,dNdxi]=lagrange_basis(coord)
  
% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%  
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
%   H4, H10, B8 and B27  
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University    

xi=coord(1); eta=coord(2);
Nv=1/4*[ (1-xi)*(1-eta);
        (1+xi)*(1-eta);
        (1+xi)*(1+eta);
        (1-xi)*(1+eta)];
dNdxi=1/4*[-(1-eta), -(1-xi);
       1-eta,    -(1+xi);
       1+eta,      1+xi;
          -(1+eta),   1-xi];