function I=trapz2(X,Y,F)
    if size(X)~=size(F)
        [X,Y]=meshgrid(X,Y);
    end
    I=trapz(Y(:,1),trapz(X(1,:),F,2));
end