%save in a structure all the variables smaller than a threshold
%04/03/2020: created
%provide max_size
 
SMALL=struct;
allvars=whos;
sizes=[allvars.bytes];
selvar=find(sizes<max_size);
for ID_var=1:length(selvar)
    if sum(strcmp(allvars(selvar(ID_var)).class,{'single','double','char','cell','struct'}))>0
        SMALL=setfield(SMALL,allvars(selvar(ID_var)).name,eval(allvars(selvar(ID_var)).name)); %#ok<SFLD>
    end
end
