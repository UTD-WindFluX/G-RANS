%% Geometrical parameters
D = [0.894  0.904]; %m
U_infty = 10;%m/s
x_T = [2    5]; %streamwise turbine locations; CHANGED TO MATCH PARABOLIC DATA
rho = 1.2; %kg/m^3
Nt = length(x_T); %number of turbines

xmin = 0; xmax = 10; %axial limits; 
rmin = 0; 
if isequal(up_BC,'Stress_free') == 1
    rmax = 2;
    x_Map = 2; Lx=5;
    r_Map = 1; Lr = 0.5;
else
    rmax = 1.2563/D(1);%radial limits; same blockage factor (0.1283)
    x_Map = 1; Lx = 5; 
    r_Map = 2; Lr = 0.5;if TSR_theor==6 Lr = 0.8; end
end
 %mapping of GLC nodes: 0---->linear mapping; 1---->algebraic mapping; 2----->another algebraic mapping

Reynolds = 103600; 
nac=[0.5*.090/D(1); 0.5*.090/D(1)+0.05];%Nacelle radius
sigma_r = 2*nac; %Changed sigma_r to test the stall zone
Nx = 90; 
esp = 2;
urf = 0.4;     %Under relaxation factor for WT forcing[0,1]
urf_tke = 0.7; %Under relaxation factor for tke production [0,1]
CDhub = 1.1;
L_diss = 0.3; %Characteristic length for Tke dissipation term
TI = 0.01; %Turbulence intensity on the inlet
% load('./BT1_DATA/tke_parabolic.mat');
tke_in = 1.5e-03;

NxGC = Nx-2;NrGC = Nr - 2;%number of GC points over x and r
invR_max = 1e10; %max. value of 1/r (the result doesn't depend on this value)
errmin = 4e-7; %max. of L2 error acceptable
maxIter = 60; %maximum number of Newton iterations
Diagnostic = 1;
Max_NuT_smooth = .8;
%% 
[Matrices] = Create_matrices(Nx,Nr,NxGC,NrGC,xmin,xmax,rmax,r_Map,x_Map,Lx,Lr);
[~,~,IWr,IWx]=normaL2matrix; %Integration matrices
Matrices.Mat_smooth_Gauss = Mat_smooth_Gauss(NuT_smooth_par);
%% Spreading matrices and initial field for Newton iterations
Ind = zeros(length(TSR),1);Ind_hub = Ind;nac_ind = Ind;
for k=1:length(TSR)
    nac_ind(k) = min(find((Matrices.RmGLC(:,1)-nac(k))<=0)); %Lowest radial index out of nacelle's radius
    Ind(k)=min(find(Matrices.XmGLC(1,:)-(x_T(k))<=0)); %index of axial induction location
    Ind_hub(k) = min(find(Matrices.XmGLC(1,:)-(x_T(k)-1)<=0)); %index for relative velocity on the nacelle
end
Ind_R(1) = min(find(Matrices.RmGLC(:,1)<=0.5)); %index of highest radial point inside Actuator Disk
Ind_R(2) = min(find(Matrices.RmGLC(:,1)<=0.55)); %index of highest radial point inside Actuator Disk

if isequal(initial_profile,'Variable')==1
    Sigma = 0.5;
    fr = @(r) 1 - exp(-0.5*r.^2./Sigma^2)./(sqrt(2*pi)*Sigma);
%     fr = @(r) r./rmax;
    start = fr(Matrices.RmGLC);
    start = reshape(start.',Nx*Nr,1);%axial velocity field
    Initial_field = [0*start; 0*start; start; zeros(NxGC*NrGC,1); tke_in*ones(Nx*Nr,1)];qc_initial = []; %Initial field (Ur Ut Ux p tke), div. free
    Inlet.Ur = Initial_field(1:Nx*Nr);
    Inlet.Ut = Initial_field(Nx*Nr+1:2*Nx*Nr);
    Inlet.Ux = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
    Inlet.ALL = Initial_field;
else
    start = ones(Nr,Nx); 
    start = reshape(start.',Nx*Nr,1);%axial velocity field
    if isequal(Turbulence_model,'TKE')==0
        Initial_field = [0*start; 0*start; start; zeros(NxGC*NrGC,1)];qc_initial = []; %Initial field (Ur Ut Ux p tke), div. free
        Inlet.Ur = Initial_field(1:Nx*Nr);
        Inlet.Ut = Initial_field(Nx*Nr+1:2*Nx*Nr);
        Inlet.Ux = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
        Inlet.ALL = Initial_field;
    else
        Initial_field = [0*start; 0*start; start; zeros(NxGC*NrGC,1); tke_in*ones(Nx*Nr,1)];qc_initial = []; %Initial field (Ur Ut Ux p tke), div. free
        Inlet.Ur = Initial_field(1:Nx*Nr);
        Inlet.Ut = Initial_field(Nx*Nr+1:2*Nx*Nr);
        Inlet.Ux = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
        Inlet.ALL = Initial_field;
    end
end
%%
airfoil = BlindTest1(Matrices.RmGLC(Ind_R:end,1));
ch = airfoil(:,2); tw = airfoil(:,3); Nfoil = airfoil(:,4);
Nrd = length(ch);
R_AD = Matrices.RmGLC(Ind_R:end,1); %radial GLC points strictly included within AD
%%
sigma_x = .3*ones(Nt,1); %st. dev. of Gaussian spreading functions along x;
sigma_t =.05*ones(Nt,1); %st. dev. of Gaussian spreading functions along theta