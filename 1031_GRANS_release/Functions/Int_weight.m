function [alpha_int] = Int_weight(xL_vett)
    NxL=length(xL_vett);
    alpha_int=zeros(size(xL_vett));

    alpha_int(1)=1/2*(xL_vett(2)-xL_vett(1));

    for nn=2:NxL-1
        alpha_int(nn)=1/2*(xL_vett(nn+1)-xL_vett(nn-1));
    end

    alpha_int(NxL)=1/2*(xL_vett(NxL)-xL_vett(NxL-1));
end
