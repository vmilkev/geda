function set_report_param( this, t_name, f_name )
    this.text_fname = t_name;
    this.figs_fname = f_name;
    if exist(this.text_fname, 'file')
        delete(this.text_fname);
    end
    if exist(this.figs_fname, 'file')
        delete(this.figs_fname);
    end
    this.f_hdl = fopen(this.text_fname,'a');
end
