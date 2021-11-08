function s=extract_filename(filename)
    filename=replace(filename,'\','/');
    slashes=regexp(filename,'/');
    dots=regexp(filename,'\.');
    s=filename((slashes(end)+1):(dots(end)-1));
end