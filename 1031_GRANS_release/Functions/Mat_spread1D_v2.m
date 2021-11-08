%08/14/2020 (v 2): truncated Guassian

function M=Mat_spread1D_v2(y,sigma,sigma_cut)
    for k=1:length(y)
       G=exp(-(y-y(k)).^2/(2*sigma.^2));   
       G(abs(y-y(k))>=sigma_cut*sigma)=0;
       I=trapz(y,G);
       M(:,k)=G/I;
    end
end