%% Numerical features
Buffer_zone                 = 0;                                            %Set 1 if EV must have a buffer zone on the last diameter; 0 otherwise
out_BC                      = 'Neumann';                                    %Choice of outlet boundary conditions
up_BC                       = 'Stress_free';                                %Set the upper boundary condition
Under_relaxation_factor_tke = 0;                                            %Set 1 if tke production must be under-relaxed during the iterations, 0 otherwise
initial_profile             = 'Constant';                                   %Starting profile for velocity
Turbulence_model            = 'ML';                                         %Turbulence model (choose between 'EV','ML','TKE');
NuT_smooth_par              = 0.0;                                          %Smoothing parameter for Eddy Viscosity (0=no smoothing, 1=max. smoothing)
Nx                          = 50; 
Nr                          = 40;

urf                         = 1;                        %Under relaxation factor about the forcing term[0,1] (keep it at 1 possibly)
urf_tke                     = 1;                        %Under relaxation factor for tke production [0,1]
CDhub                       = 1.1;                      %Drag coefficient of the nacelle
L_diss                      = 0.3;                      %Characteristic length for Tke dissipation term
TI                          = 0.01;                     %Turbulence intensity on the inlet
tke_in                      = 1.5*TI^2;                 %Turbulent kinetic energy at the inlet

NxGC                        = Nx-2;
NrGC                        = Nr - 2;                   %number of GC points over x and r
invR_max                    = 1e10;                     %max. value of 1/r (the result doesn't depend on this value)
errmin                      = 4e-7;                     %max. of acceptable L2 error 
maxIter                     = 60;                       %maximum number of Newton iterations
Diagnostic                  = 0;                        %Set 1 if velocity must be visualized during the iterations
Variation                   = 1;                        %Set 1 if the increment of velocity must be visualized during the iterations;
Max_NuT_smooth              = .8;                       %Maximim value of EV smoothing parameter
%% Geometrical parameters
D           = 0.894;                    %m, turbine diameted
U_infty     = 10;                       %m/s, free stream velocity
x_T         = 2;                        %streamwise turbine locations;
Nt          = length(x_T);              %number of turbines

xmin        = 0; 
xmax        = 7;                        %axial limits;
rmin        = 0; 
TSR         = 5;
Reynolds    = U_infty*D/1.45e-05;       %Reynolds number (U_infty D/nu)
nac         = 0.5*.090/D;               %Nacelle radius
sigma_r     = 2*nac;                    %Spreading standard deviation for the nacelle drag

if isequal(up_BC,'Stress_free') == 1
    rmax    = 1.5;
    x_Map   = 2; 
    Lx      = 5;
    r_Map   = 1; 
    Lr      = 0.5;
else
    rmax    = 1.2563/D;                 %radial limits; same blockage factor (0.1283)
    x_Map   = 2; 
    Lx      = 5;          %
    r_Map   = 2; 
    Lr      = 0.5;
end
%mapping of GLC nodes: 0---->linear mapping; 1---->algebraic mapping;
%2----->algebraic mapping clustering nodes in the middle of the domain
%% Create differentiation matrices
[Matrices]                  = Create_matrices(Nx,Nr,NxGC,NrGC,xmin,xmax,rmax,r_Map,x_Map,Lx,Lr,9);   %Create data structure with geometric info
nac_ind                     = find((Matrices.RmGLC(:,1)-nac)<=0,1);                                  %Lowest radial index out of nacelle's radius
[~,~,IWr,IWx]               = normaL2matrix;                                                        %Integration matrices
Matrices.Mat_smooth_Gauss   = Mat_smooth_Gauss(NuT_smooth_par);                     %Include matrix with smoothing parameter
%% Spreading matrices and initial field for Newton iterations
Ind         = zeros(length(TSR),1);
Ind_hub     = Ind;                                                                  %Initialization of indeces for turbine location

for k=1:length(TSR)
    Ind(k)      = find(Matrices.XmGLC(1,:)-(x_T(k))<=0,1);                          %index of axial induction location
    Ind_hub(k)  = find(Matrices.XmGLC(1,:)-(x_T(k)-1)<=0,1);                        %index for relative velocity on the nacelle
end
Ind_R       = find(Matrices.RmGLC(:,1)<=0.5,1);                                     %index of highest radial point inside Actuator Disk

if isequal(initial_profile,'Variable')                                              %initial_profile = Variable imposes a Gaussian profile at the inlet
    Sigma               = 0.5;
    fr                  = @(r) 1 - exp(-0.5*r.^2./Sigma^2)./(sqrt(2*pi)*Sigma);     %Function handle for Gaussian

    start               = fr(Matrices.RmGLC);                                       %Evaluate Gaussian on radial nodes
    start               = reshape(start.',Nx*Nr,1);                                 %axial velocity field
    Initial_field       = [0*start; 0*start; start; zeros(NxGC*NrGC,1); tke_in*ones(Nx*Nr,1)];
    qc_initial          = []; %Initial field (Ur Ut Ux p tke), div. free
    Inlet.Ur            = Initial_field(1:Nx*Nr);                                              %Initialization data structure for initial field
    Inlet.Ut            = Initial_field(Nx*Nr+1:2*Nx*Nr);
    Inlet.Ux            = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
    Inlet.ALL           = Initial_field;
else
    start = ones(Nr,Nx);                                                            %Constant velocity field (if chosen by the user)
    start = reshape(start.',Nx*Nr,1);                                               %axial velocity field
    if isequal(Turbulence_model,'TKE')
        Initial_field = [0*start; 0*start; start; zeros(NxGC*NrGC,1); tke_in*ones(Nx*Nr,1)];qc_initial = []; %Initial field (Ur Ut Ux p tke), div. free
        Inlet.Ur = Initial_field(1:Nx*Nr);
        Inlet.Ut = Initial_field(Nx*Nr+1:2*Nx*Nr);
        Inlet.Ux = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
        Inlet.ALL = Initial_field;
    else
        Initial_field = [0*start; 0*start; start; zeros(NxGC*NrGC,1)];qc_initial = []; %Initial field (Ur Ut Ux p tke), div. free
        Inlet.Ur = Initial_field(1:Nx*Nr);
        Inlet.Ut = Initial_field(Nx*Nr+1:2*Nx*Nr);
        Inlet.Ux = Initial_field(2*Nx*Nr+1:3*Nx*Nr);
        Inlet.ALL = Initial_field;
    end
end
%% Arifoil geometry
airfoil     = BlindTest1(Matrices.RmGLC(Ind_R:end,1));                                  %Recover geometry of the blade
ch          = airfoil(:,2); 
tw          = airfoil(:,3); 
Nfoil       = airfoil(:,4);                                                             %ch = Chord on spanwise; tw=twist angle on spanwise; Nfoil = number of airfoils
Nrd         = length(ch);
R_AD        = Matrices.RmGLC(Ind_R:end,1);                                              %radial GLC points strictly included within AD
%% Gaussian spreading forcing
sigma_x = .3*ones(Nt,1);                                                            %st. dev. of Gaussian spreading functions along x;
sigma_t = .1*ones(Nt,1);                                                            %st. dev. of Gaussian spreading functions along theta