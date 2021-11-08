function [JA, F] = JacobianSponge_fast_forced(qc,nut,DnutDUr,DnutDUth,DnutDUx,Forcing)
Globals

D1_X=Matrices.D1_X;
D2_X=Matrices.D2_X;
D1_r=Matrices.D1_r;
D2_r=Matrices.D2_r;
D1GC_X=Matrices.D1GC_X;
D1GC_r=Matrices.D1GC_r;
MfrL=Matrices.MfrL;
MtoL=Matrices.MtoL;
RmGLC=Matrices.RmGLC;
RmGC=Matrices.RmGC;


%%

rGLCc=reshape(RmGLC',Nx*Nr,1);
rGCc=reshape(RmGC',NxGC*NrGC,1);

Urc=qc(1:Nx*Nr);
Ur  = diag(Urc);
DUrDX=diag(D1_X*Urc);
DUrDr=diag(D1_r*Urc);

Utc=qc(Nx*Nr+1:2*Nx*Nr);
DUtDX=diag(D1_X*Utc);
DUtDr=diag(D1_r*Utc);

Uxc=qc(2*Nx*Nr+1:3*Nx*Nr);
Ux  = diag(Uxc);
DUxDX=diag(D1_X*Uxc);
DUxDr=diag(D1_r*Uxc);

pc = qc(3*Nx*Nr+1:end);

%
invrGLCc=1./rGLCc; invrGLCc(end-Nx+1:end)=invR_max; %the result doesn't depend on this value
invrGLC=diag(invrGLCc);

Om=diag(Utc.*invrGLCc); Om(Nx*Nr-Nx+1:end,:)=0;

nuRe=diag(1/Reynolds+nut);
DnutDx=diag(Matrices.D1_X*nut); 
DnutDr=diag(Matrices.D1_r*nut); 

Gamma=Ur*D1_r+Ux*D1_X;
Delta=invrGLC*D1_r+D2_r+D2_X;
Delta2=Delta-diag(invrGLCc.^2);


% JA
JA11=DUrDr ;
JA12=-Om;
JA13=DUrDX;
JA14=zeros(Nx*Nr,NxGC*NrGC); 

JA21=DUtDr;
JA22=invrGLC*Ur;
JA23=DUtDX;
JA24=zeros(Nx*Nr,NxGC*NrGC);

JA31=DUxDr;
JA32=zeros(Nr*Nx);
JA33=DUxDX;
JA34=zeros(Nx*Nr,NxGC*NrGC);

JA4=zeros(NxGC*NrGC,Nx*Nr*3+NxGC*NrGC);



% FA
FA11=Gamma - nuRe*Delta2  -DnutDx*D1_X- 2*DnutDr*D1_r;
FA12=-Om;
FA13= -DnutDx*D1_r;
FA14=MtoL*D1GC_r;

FA21=Om;
FA22=Gamma-nuRe*Delta2  -DnutDx*D1_X-DnutDr*(-invrGLC + D1_r);
FA23=zeros(Nx*Nr);
FA24=zeros(Nx*Nr,NxGC*NrGC);

FA31=-DnutDr*D1_X;
FA32=zeros(Nx*Nr);
FA33=(Gamma-nuRe*Delta) -2*DnutDx*D1_X- DnutDr*D1_r;
FA34=MtoL*D1GC_X;

FA41=diag(1./rGCc)*MfrL+MfrL*D1_r;
FA42=zeros(NxGC*NrGC,Nx*Nr);
FA43=MfrL*D1_X;
FA44=zeros(NxGC*NrGC);



J_nl =  [JA11, JA12, JA13, JA14; JA21, JA22, JA23, JA24; JA31, JA32, JA33, JA34; JA4];
FA = [FA11, FA12, FA13, FA14; FA21, FA22, FA23, FA24; FA31, FA32, FA33, FA34; FA41, FA42, FA43, FA44];

switch Turbulence_model
    case 'ML'
        DR=(diag(Delta2*Urc) + 2*DUrDr*D1_r + (DUrDX + DUxDr)*D1_X);
        DT=(diag(Delta2*Utc)+DUtDX*D1_X+diag((D1_r-invrGLC)*Utc)*D1_r);
        DX=(diag(Delta*Uxc)+2*DUxDX*D1_X+(DUrDX + DUxDr)*D1_r);

        %% Jacobian blocks that must be joined to the others due to the ML closure
        JA11= -DR*DnutDUr;
        JA12= -DR*DnutDUth;
        JA13= -DR*DnutDUx;
        JA14= zeros(Nx*Nr,NxGC*NrGC);

        JA21= -DT*DnutDUr;
        JA22= -DT*DnutDUth;
        JA23= -DT*DnutDUx;
        JA24= zeros(Nx*Nr,NxGC*NrGC);

        JA31= -DX*DnutDUr;
        JA32= -DX*DnutDUth;
        JA33= -DX*DnutDUx;
        JA34= zeros(Nx*Nr,NxGC*NrGC);

        JA4=zeros(NxGC*NrGC,Nx*Nr*3+NxGC*NrGC);

        JA_adj =  [JA11, JA12, JA13, JA14; JA21, JA22, JA23, JA24; JA31, JA32, JA33, JA34; JA4];

    case 'EV'
        JA_adj=0;
end

F=FA*qc-Forcing;
A=FA+J_nl+JA_adj;

%% boundary conditions

col_Ur=0*Nx*Nr+1:Nx*Nr*1;
col_Ut=1*Nx*Nr+1:Nx*Nr*2;
col_Ux=2*Nx*Nr+1:Nx*Nr*3;
col_p= 3*Nx*Nr+1:Nx*Nr*3+NxGC*NrGC;

%% top
%stress free
row_TOP=1:Nx;
pGLC=MtoL*pc;

F(row_TOP) = 2*diag(nuRe(row_TOP,row_TOP)*DUrDr(row_TOP,row_TOP))-pGLC(row_TOP); %2muuu durdr-p=0 ---> (2/Re) durdr-p=0

A(row_TOP,col_Ur)=2*nuRe(row_TOP,row_TOP)*D1_r(row_TOP,:);
A(row_TOP,col_Ut)=0;
A(row_TOP,col_Ux)=0;
A(row_TOP,col_p)=-MtoL(row_TOP,:);

% 
F(row_TOP+Nr*Nx)=(D1_r(row_TOP,:)+invrGLC(row_TOP,:))*Utc; % -1/r(uth)+duthdr=0

A(row_TOP+Nr*Nx,col_Ur)=0;
A(row_TOP+Nr*Nx,col_Ut)=(D1_r(row_TOP,:)+invrGLC(row_TOP,:));
A(row_TOP+Nr*Nx,col_Ux)=0;
A(row_TOP+Nr*Nx,col_p)=0;

% 
F(row_TOP+2*Nr*Nx)=diag(DUrDX(row_TOP,row_TOP)+DUxDr(row_TOP,row_TOP)); % durdx+duxdr=0

A(row_TOP+2*Nr*Nx,col_Ur)=D1_X(row_TOP,:);
A(row_TOP+2*Nr*Nx,col_Ut)=0;
A(row_TOP+2*Nr*Nx,col_Ux)=D1_r(row_TOP,:);
A(row_TOP+2*Nr*Nx,col_p)=0;



 %% outlet
    row_OUTLET=1:Nx:Nx*Nr;
    switch out_BC
        case 'Neumann'
            %homogeneous neumann

            F(row_OUTLET)=diag(DUrDX(row_OUTLET,row_OUTLET)); 

            A(row_OUTLET,col_Ur)=D1_X(row_OUTLET,:);
            A(row_OUTLET,col_Ut)=0;
            A(row_OUTLET,col_Ux)=0;
            A(row_OUTLET,col_p)=0;
            %
            F(row_OUTLET+Nx*Nr)=diag(DUtDX(row_OUTLET,row_OUTLET)); 

            A(row_OUTLET+Nx*Nr,col_Ur)=0;
            A(row_OUTLET+Nx*Nr,col_Ut)=D1_X(row_OUTLET,:);
            A(row_OUTLET+Nx*Nr,col_Ux)=0;
            A(row_OUTLET+Nx*Nr,col_p)=0;

            %
            F(row_OUTLET+2*Nx*Nr)=diag(DUxDX(row_OUTLET,row_OUTLET)); 

            A(row_OUTLET+2*Nx*Nr,col_Ur)=0;
            A(row_OUTLET+2*Nx*Nr,col_Ut)=0;
            A(row_OUTLET+2*Nx*Nr,col_Ux)=D1_X(row_OUTLET,:);
            A(row_OUTLET+2*Nx*Nr,col_p)=0;
        case 'non-homogeneous_Neumann'
            %non-homogeneous neumann
            I=eye(Nx*Nr);

            F(row_OUTLET)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*diag(DUrDX); 

            A(row_OUTLET,col_Ur)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*D1_X;
            A(row_OUTLET,col_Ut)=0;
            A(row_OUTLET,col_Ux)=0;
            A(row_OUTLET,col_p)=0;
            %
            F(row_OUTLET+Nx*Nr)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*diag(DUtDX); 

            A(row_OUTLET+Nx*Nr,col_Ur)=0;
            A(row_OUTLET+Nx*Nr,col_Ut)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*D1_X;
            A(row_OUTLET+Nx*Nr,col_Ux)=0;
            A(row_OUTLET+Nx*Nr,col_p)=0;

            %
            F(row_OUTLET+2*Nx*Nr)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*diag(DUxDX); 

            A(row_OUTLET+2*Nx*Nr,col_Ur)=0;
            A(row_OUTLET+2*Nx*Nr,col_Ut)=0;
            A(row_OUTLET+2*Nx*Nr,col_Ux)=(I(row_OUTLET,:)-I(row_OUTLET+1,:))*D1_X;
            A(row_OUTLET+2*Nx*Nr,col_p)=0;
    end

    %
    row_OUTLET_GC=1:NxGC:NxGC*NrGC;
    I_GC=eye(NxGC*NrGC);
    
    F(row_OUTLET_GC+3*Nx*Nr)=pc(row_OUTLET_GC); 
    
    A(row_OUTLET_GC+3*Nx*Nr,col_Ur)=0;
    A(row_OUTLET_GC+3*Nx*Nr,col_Ut)=0;
    A(row_OUTLET_GC+3*Nx*Nr,col_Ux)=0;
    A(row_OUTLET_GC+3*Nx*Nr,col_p)=I_GC(row_OUTLET_GC,:);


%% axis 
row_AXIS=Nr*Nx-Nx+1:Nr*Nx;
% according to Batchelor&Gill 1962  (e come Meyer&Powell 1992 e Olendraru&Sellier 2002)

F(row_AXIS)=Urc(row_AXIS); 

I=eye(Nr*Nx);

A(row_AXIS,col_Ur)=I(row_AXIS,:);
A(row_AXIS,col_Ut)=0;
A(row_AXIS,col_Ux)=0;
A(row_AXIS,col_p)=0;

F(row_AXIS+Nx*Nr)=Utc(row_AXIS); 

A(row_AXIS+Nx*Nr,col_Ur)=0;
A(row_AXIS+Nx*Nr,col_Ut)=I(row_AXIS,:);
A(row_AXIS+Nx*Nr,col_Ux)=0;
A(row_AXIS+Nx*Nr,col_p)=0;

F(row_AXIS+2*Nx*Nr)=diag(DUxDr(row_AXIS,row_AXIS)); 

A(row_AXIS+2*Nx*Nr,col_Ur)=0;
A(row_AXIS+2*Nx*Nr,col_Ut)=0;
A(row_AXIS+2*Nx*Nr,col_Ux)=D1_r(row_AXIS,:);
A(row_AXIS+2*Nx*Nr,col_p)=0;


%% inlet
%dirichlet

row_INLET=Nx:Nx:Nr*Nx;

F(row_INLET)=Urc(row_INLET)-Inlet.ALL(row_INLET); 

A(row_INLET,col_Ur)=I(row_INLET,:);
A(row_INLET,col_Ut)=0;
A(row_INLET,col_Ux)=0;
A(row_INLET,col_p)=0;

F(row_INLET+Nx*Nr)=Utc(row_INLET)-Inlet.ALL(row_INLET+Nx*Nr); 

A(row_INLET+Nx*Nr,col_Ur)=0;
A(row_INLET+Nx*Nr,col_Ut)=I(row_INLET,:);
A(row_INLET+Nx*Nr,col_Ux)=0;
A(row_INLET+Nx*Nr,col_p)=0;

F(row_INLET+2*Nx*Nr)=Uxc(row_INLET)-Inlet.ALL(row_INLET+2*Nx*Nr); 

A(row_INLET+2*Nx*Nr,col_Ur)=0;
A(row_INLET+2*Nx*Nr,col_Ut)=0;
A(row_INLET+2*Nx*Nr,col_Ux)=I(row_INLET,:);
A(row_INLET+2*Nx*Nr,col_p)=0;


JA=A;
