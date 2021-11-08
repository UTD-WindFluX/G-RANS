%Newton solver of the RANS field

function [RANS]=AxiSymDirect_EV(NuT_coeff,Forcing)
    Globals
    NuT=NuT_fun(NuT_coeff);
    NuTc=reshape(NuT',Nx*Nr,1);
    if isempty(qc_initial) %it means that this is the first iteration
        qc0=Inlet.ALL;
    else
        qc0=qc_initial;
    end

    SquartoL2=1e5*errmin;%norm initialization
    nIte=0;
    UR_forcing=0.1;
    while (SquartoL2>errmin || UR_forcing<1)
        if nIte==0;qc=qc0;end
        
        [JA,F] = JacobianSponge_fast_forced(qc,NuTc,0,0,0,Forcing*UR_forcing);%jacobain and residuals calculator

        deltaq=-JA\F;%unknowns variation
        qc=qc+deltaq;%new guess
         
        [SquartoL2]=normaL2(deltaq);%norm of the unknown vector variation
        display(['SqL2=',num2str(SquartoL2)]);    

        nIte=nIte+1;
        UR_forcing=min(nIte/3,1);
        Error(nIte) = SquartoL2;
        if SquartoL2>10 || nIte>maxIter
            disp('No convergence');
            go_on=false;
            RANS=0;
            return
        end   
    end
   
    %% Rearrangment of the results
    Ur=reshape(qc(1:Nx*Nr),Nx,Nr)'; Ut=reshape(qc(Nx*Nr+1:2*Nx*Nr),Nx,Nr)';
    Ux=reshape(qc(2*Nx*Nr+1:3*Nx*Nr),Nx,Nr)'; p=reshape(Matrices.MtoL*qc(3*Nx*Nr+1:3*Nx*Nr+NxGC*NrGC),Nx,Nr)';
    Fx=reshape(Forcing(2*Nx*Nr+1:3*Nx*Nr)*UR_forcing,Nx,Nr)';   
    
    %% output
    RANS.Ur=Ur;RANS.Ut=Ut;RANS.Ux=Ux;RANS.p=p;
    RANS.L2Err=Error;RANS.NuT=NuT;           
    RANS.Fx=Fx;
