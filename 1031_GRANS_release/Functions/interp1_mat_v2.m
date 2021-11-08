%07/31/2020 (v 2): added nans
%05/13/2021: finalized
function M=interp1_mat_v2(source_grid, target_grid)
    M=zeros(length(target_grid),length(source_grid));
    for i=1:length(target_grid)
        if source_grid(end)>source_grid(1) 
            j2=find(source_grid>=target_grid(i),1,'first');
            j1=find(source_grid<=target_grid(i),1,'last');
        else
            j2=find(source_grid>=target_grid(i),1,'last');
            j1=find(source_grid<=target_grid(i),1,'first');
        end
        x1=source_grid(j1);x2=source_grid(j2);x=target_grid(i);
        if isempty(j1) || isempty(j2)
            M(i,:)=0;
        else
            if j1==j2
                M(i,j1)=1;
            else
                M(i,j1)=1-(x-x1)/(x2-x1+eps);
                M(i,j2)=(x-x1)/(x2-x1+eps);
            end
        end
    end
    M=sparse(M);
end