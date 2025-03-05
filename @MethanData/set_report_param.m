function set_report_param( this, t_name )
    this.text_fname = t_name;
    if exist(this.text_fname, 'file') % though, this is unlikely event, but is possible
        r = randi(1000,1);
        disp("WARNING: There is the log file with the same name already exists in the working folder.");
        this.text_fname = strcat("log_geda_",this.date_now,"_",num2str(r),".txt");
        disp( strcat("         Therefore the log file name will be modified to: ", this.text_fname) );
        %delete(this.text_fname); will not delete !
    end
    this.f_hdl = fopen(this.text_fname,'a');
    if this.f_hdl == -1
        warning("Cannot open the program's log file!");
        return;
    end
    this.make_report( "dat", strcat("GEDA LOG, starts: ",string(datetime()) ), []);
    this.make_report("dat", " ", []);
end
