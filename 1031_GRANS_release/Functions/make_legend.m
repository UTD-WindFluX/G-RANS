function l=make_legend(x,name,unit,type)
    switch type
        case 'center'
            for i=1:length(x)
                l{i}=[name,' $ = ',num2str(x(i)),'$',unit];
            end
        case 'bins'
            for i=1:length(x)-1
                l{i}=[name,' $\in [',num2str(x(i)),', ',num2str(x(i+1)),')$',unit];
            end
    end
end