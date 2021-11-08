function save_fig(name)
    f=gcf;
    ax=gca;
    
    c=clock;
    
    day=c(3);
    month=c(2);
    folder=['C:\Users\User\Desktop\PhD\Analysis\',num2str(month),'_',num2str(day),'\'];
    if exist(folder,'file')~=7
        mkdir(folder)
    end
    set(f,'InvertHardCopy','off');
    saveas(f,[folder,name,'.png']);
    savefig(f,[folder,name,'.fig']);
  
    
end