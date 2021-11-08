function f=trunc_Gauss(x,sigma,sigmas_cut)
    f=exp(-x.^2/(2*sigma^2));
    f(abs(x)>sigmas_cut*sigma)=0;
    f=f/trapz(x,f);
end