function F=fade_2D(X,Y,Lx,Ly)
    fx=ones(size(X));
    fy=ones(size(Y));
    
    sx=sign(X(1,2)-X(1,1));
    sy=sign(Y(2,1)-Y(1,1));
    
    sel=(X-X(1,1))*sx<Lx;
    fx(sel)=3*(X(sel)-X(1,1)).^2/Lx^2-2*(X(sel)-X(1,1)).^3*sx/Lx^3;
    sel=(X(1,end)-X)*sx<Lx;
    fx(sel)=3*(X(1,end)-X(sel)).^2/Lx^2-2*(X(1,end)-X(sel)).^3*sx/Lx^3;
    
    sel=(Y-Y(1,1))*sy<Ly;
    fy(sel)=3*(Y(sel)-Y(1,1)).^2/Ly^2-2*(Y(sel)-Y(1,1)).^3*sy/Ly^3;
    sel=(Y(end,1)-Y)*sy<Ly;
    fy(sel)=3*(Y(end,1)-Y(sel)).^2/Ly^2-2*(Y(end,1)-Y(sel)).^3*sy/Ly^3;
    
    F=fx.*fy;
end