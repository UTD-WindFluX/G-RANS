function [Squarto]=normaL2(deltaU)
Globals

x1GLC=Matrices.x1GLC;
rGLC=Matrices.RmGLC(:,1);
y1GLC=Matrices.y1GLC;


IWxLxW=diag(INTweights(Nx,2)); 
IWxMLxW=diag(IWxLxW./x1GLC); 

IWy=diag(INTweights(Nr,2));
IWyM=diag(IWy.*rGLC./y1GLC);

IWMLxW=kron(IWyM,IWxMLxW);

Z=0*IWMLxW;

PE=[IWMLxW,Z,Z;Z,IWMLxW,Z;Z,Z,IWMLxW];
Squarto=sqrt(deltaU(1:3*Nx*Nr)'*PE*deltaU(1:3*Nx*Nr));
