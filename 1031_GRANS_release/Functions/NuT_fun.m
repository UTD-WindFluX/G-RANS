%This function calculates the EV as a function of x for given coefficent(s)
%1)if the coefficent is a scalar, the EV is constant and equal to the coefficient
%2)if the coefficients are a 2x1 array [A B], the EV is equal to A+Bx, and the x is zero at the turbine's hub 
function NuT=NuT_fun(NuT_coeff)
    Globals
   
    switch DoFs
        case 1
            NuT=ones(Nr,Nx)*NuT_coeff(1);
        case 2
            NuT=(XmGLC+xmin)*NuT_coeff(2)+NuT_coeff(1);
    end
    
end