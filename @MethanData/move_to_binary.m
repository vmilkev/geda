function move_to_binary( this, data, fname )
    fileID = fopen(fname,'w');
    if fileID == -1
        warning("fopen cannot open the file!");
        return;
    end
    fwrite(fileID, size(data), 'int32');
    fwrite(fileID,data,'double');
    fclose(fileID);
end