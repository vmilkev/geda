function fclean( this, files_list )
    for i = 1:size(files_list,1)
        delete(files_list{i,1});
    end
end