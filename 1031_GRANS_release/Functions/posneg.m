function map=posneg(m)
    narginchk(0,1);  

    if nargin == 1
        validateattributes(m,{'numeric'},{'numel',1});
    end
    if nargin < 1, m = size(get(gcf,'colormap'),1); 
    end
    c=[0 0 0.5;0 0.5 1; 1 1 1; 1 0 0; 0.5 0 0];
    pp=1:(m-1)/(size(c,1)-1):m;
    r=interp1(pp,c(:,1),1:m);
    g=interp1(pp,c(:,2),1:m);
    b=interp1(pp,c(:,3),1:m);

    map=[r' g' b']/255;
    map = map/max(map(:));
end