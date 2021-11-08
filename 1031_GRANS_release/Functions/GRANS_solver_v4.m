%GRANS single turbine solver
%06/04/2021 (v 2): load externalized
%09/19/2021 (v 3): loads directly from inverse depth-average
%09/20/2021: finalized
%11/01/2021 (v 4): adapted for sharing, finalized

Globals

%% Inputs
Turbulence_model='EV';
Reynolds=10^7;
DoFs=1;%degrees of freedom of turbulence model
BL_approx=0;%use boundary laywer approximation
spr_load=0.5;%SD of the Gaussian spreading of loads
NuT_smooth_par=0.1;%sigma of the gaussian smoothing fucntion for EV (used only with ML model)
Max_NuT_smooth=0.5;%limit of the smoothing parameter
Diagnostic=0;%display detailed plots during Newton loop
invR_max=10^10;%maximum value of invrGLC
errmin=10^-8;%maximum allowable error for the Newton method
maxIter=30; %maximum allowable number of Newton iterations
out_BC='non-homogeneous_Neumann';

%% Initiliazation

%load input file
Ct=xlsread(source_input,'GRANS input','B1');
NuT=xlsread(source_input,'GRANS input','B2');
xmin=xlsread(source_input,'GRANS input','B3');
xmax=xlsread(source_input,'GRANS input','B4');
rmax=xlsread(source_input,'GRANS input','B5');
load_r=xlsread(source_input,'GRANS input','D2:D500');
load_phi=xlsread(source_input,'GRANS input','E2:E500');

% GLC domain
L_domain=xmax-xmin;%xmax-xmin
x_Map=1;%mapping along x: 0--> linear; 1--> algebraic
r_Map=1;%mapping along r: 0--> linear; 1--> algebraic
Lx=0.5*L_domain;%algebraic x mapping parameter
Lr=0.5*rmax;%algebraic r-mapping parameter
Nx=round(5*L_domain); %Number of GLC points along x
Nr=max(round(15*rmax),20); %Number of GLC points along r
NxGC=Nx-2; %number of points in the GC grid along x axis
NrGC=Nr-2; %number of points in the GC grid along r axis
Create_matrices%creates matrices for derivation and grid-change

%zeroing
go_on=true;
Inlet.ALL=[zeros(2*Nx*Nr,1);ones(Nr*Nx,1);zeros(NxGC*NrGC,1)];
qc_initial=[]; 

%forcing
c=Ct/(16*trapz(load_r,load_phi.*load_r));
s=exp(-(xmin+XmGLC(1,:)).^2/(2*spr_load^2));
s=s/abs(trapz(xmin+XmGLC(1,:),s));
phi=interp1(load_r,load_phi,RmGLC(:,1));
phi(isnan(phi))=0;
Fx=-c*s(:)'.*phi(:);
Forcing=[zeros(Nx*Nr*2,1);vert(Fx') ;zeros(NxGC*NrGC,1)];

%% Main
AXI=AxiSymDirect_EV(NuT,Forcing);%RANS solver
AXI.x=XmGLC(1,:)+xmin;
AXI.r=RmGLC(:,1);
AXI.Ct=Ct;
AXI.NuT=NuT;
AXI.xmin=xmin;
AXI.xmax=xmax;
AXI.rmax=rmax;
AXI.load_r=load_r;
AXI.load_phi=load_phi;
AXI.source_input=source_input;

%% Save
max_size=1000;find_small;
save(['./Results/',extract_filename(source_input),'.mat'],'AXI');

