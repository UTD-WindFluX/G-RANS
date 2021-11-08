function M=Mat_spread1D(y,sigma)
    for k=1:length(y)
       G=exp(-(y-y(k)).^2/(2*sigma.^2));     
       I=trapz(y,G);
       M(:,k)=G/I;
    end
end