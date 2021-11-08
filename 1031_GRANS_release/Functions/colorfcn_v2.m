function col=colorfcn_v2(x,xmin,xmax,color_map)
    cols=feval(color_map,100);
    for i=1:3
        col(i)=interp1(0:1/(length(cols(:,1))-1):1,cols(:,i),(x-xmin)/(xmax-xmin));
    end
end