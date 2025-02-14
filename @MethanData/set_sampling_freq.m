function set_sampling_freq( this, freq )
    if freq < 1.0
        this.make_report( "dat", "Warning: The provided sampling frequency cannot be 0 or negative!", [] );
    end
    this.sampl_res_const = freq;
end