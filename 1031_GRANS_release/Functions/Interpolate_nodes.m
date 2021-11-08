function [M]=Interpolate_nodes(N,GLC_points,GC_points,ind)
%%%%% ind=0----> from GLC
%%%%% ind=1----> to GLC
    PhiGLC = ChebyshevInterpolants(N,2,GLC_points);
    PhiGC = ChebyshevInterpolants(N,2,GC_points);

    if ind==0
        M=PhiGC/PhiGLC;
    elseif ind==1
        M=PhiGLC/PhiGC;
    end
end