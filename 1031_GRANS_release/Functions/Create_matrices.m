    %creates the derivation and grid-change matrices
    [xGLC,D1x]=nodiD1chebGLC(Nx-1);
    [yGLC,D1y]=nodiD1chebGLC(Nr-1);  

    [xGC,D1xGC]=nodichebGC(NxGC);
    [yGC,D1yGC]=nodichebGC(NrGC);
    
    switch r_Map
        case 0 %linear
            r_GLC=rmax/2*(1+yGLC);
            r_GC=rmax/2*(1+yGC);
            MD1r=2/rmax*D1y;
            MD1rGC=2/rmax*D1yGC;
            D2y=D1y*D1y;
            y1GLC=(2/rmax)*ones(size(r_GLC));
%             y1GC=(2/rmax)*ones(size(rGC));
            MD2r=(2/rmax)^2*D2y;

        case 1 %algebraic
            alfa=(2*Lr+rmax)/rmax;
            r_GLC=Lr*(1+yGLC)./(alfa-yGLC);
            r_GC=Lr*(1+yGC)./(alfa-yGC);
            
            y1GLC=Lr*(alfa+1)./(r_GLC+Lr).^2;
            y2GLC=-2*Lr*(alfa+1)./(r_GLC+Lr).^3;
            y1GC=Lr*(alfa+1)./(r_GC+Lr).^2;

            MD1r=diag(y1GLC)*D1y;
            MD1rGC=diag(y1GC)*D1yGC;
            D2y=D1y*D1y;
            MD2r=(diag(y1GLC)).^2*D2y+diag(y2GLC)*D1y;

    end

    switch x_Map
        case 0 %linear
            x_GLC=(xmax-xmin)/2*(1+xGLC);
            x_GC=(xmax-xmin)/2*(1+xGC);

            MD1X=2/(xmax-xmin)*D1x;
            MD1XGC=2/(xmax-xmin)*D1xGC;
            D2x=D1x*D1x;
            x1GLC=(2/(xmax-xmin))*ones(size(x_GLC));
%             x1GC=(2/xmax)*ones(size(xGC));
            MD2X=(2/(xmax-xmin))^2*D2x;

        case 1 %algebraic
            alfa=(2*Lx+(xmax-xmin))/(xmax-xmin);
            x_GLC=Lx*(1+xGLC)./(alfa-xGLC);
            x_GC=Lx*(1+xGC)./(alfa-xGC);

            x1GLC=Lx*(alfa+1)./(x_GLC+Lx).^2;
            x2GLC=-2*Lx*(alfa+1)./(x_GLC+Lx).^3;
            x1GC=Lx*(alfa+1)./(x_GC+Lx).^2;

            MD1X=diag(x1GLC)*D1x;
            MD1XGC=diag(x1GC)*D1xGC;
            D2x=D1x*D1x;
            MD2X=(diag(x1GLC)).^2*D2x+diag(x2GLC)*D1x;
    end

    %interpol
    % from Gauss lobatto (velocity) to Gauss (pressure)
    [MfrLx]=Interpolate_nodes(Nx,xGLC,xGC,0);
    [MfrLy]=Interpolate_nodes(Nr,yGLC,yGC,0);
    % from Gauss (pressure) to Gauss lobatto (velocity)
    [MtoLx]=Interpolate_nodes(Nx,xGLC,xGC,1);
    [MtoLy]=Interpolate_nodes(Nr,yGLC,yGC,1);

    %% 2D

    %GLC
    IxGLC=eye(Nx);
    IyGLC=eye(Nr);

    D1bi_X=kron(IyGLC,MD1X);
    D2bi_X=kron(IyGLC,MD2X);
    D1bi_r=kron(MD1r,IxGLC);
    D2bi_r=kron(MD2r,IxGLC);    
    [XmGLC,RmGLC]=meshgrid(x_GLC,r_GLC);

    %GC
     IxGC=eye(NxGC);
     IyGC=eye(NrGC);

     D1biGC_X=kron(IyGC,MD1XGC);
     D1biGC_r=kron(MD1rGC,IxGC);
     [XmGC,RmGC]=meshgrid(x_GC,r_GC); 

    %grid-change matrices
    MfrL=kron(MfrLy,MfrLx);
    MtoL=kron(MtoLy,MtoLx);
    
%    %% plot
%    figure;
% plot(XmGLC+xmin,RmGLC,'.r');hold on;plot(XmGC+xmin,RmGC,'.b');title('Spatial discretization');%legend('Momentum nodes','Pressure-Continuity nodes');
% plot_turbine_fcn(0,0,0,180,1);xlabel('x/D');ylabel('r/D');
% 
% axis equal; axis tight;
    %% output
    Matrices.D1_X=D1bi_X;
    Matrices.D2_X=D2bi_X;
    Matrices.D1_r=D1bi_r;
    Matrices.D2_r=D2bi_r;
    Matrices.D1GC_X=D1biGC_X;
    Matrices.D1GC_r=D1biGC_r;
    Matrices.MfrL=MfrL;
    Matrices.MtoL=MtoL;
    Matrices.XmGLC=XmGLC;
    Matrices.RmGLC=RmGLC;
    Matrices.XmGC=XmGC;
    Matrices.RmGC=RmGC;
    Matrices.x1GLC=x1GLC;
    Matrices.y1GLC=y1GLC;
%     Matrices.M_smooth=Mat_smooth_Gauss(NuT_smooth_par);%smoothing matrix
