function set_report_param( this, t_name )
    this.text_fname = t_name;
    if exist(this.text_fname, 'file')
        delete(this.text_fname);
    end
    this.f_hdl = fopen(this.text_fname,'a');
    if this.f_hdl == -1
        warning("Cannot open the program's log file!");
        return;
    end
    this.make_report( "dat", strcat("GEDA LOG, starts: ",string(datetime()) ), []);
    this.make_report("dat", " ", []);
end
