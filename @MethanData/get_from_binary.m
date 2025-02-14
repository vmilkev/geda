function data = get_from_binary( this, fname )
    fileID = fopen(fname,'r');
    if fileID == -1
        warning("fopen cannot open the file!");
        return;
    end
    size_data = fread(fileID,[1 2],'int32');
    data = fread(fileID, size_data, 'double');
    fclose(fileID);
end