function read_lely_data( this )

[lely_map, lely_num_map, robots] = this.importfile_xlsx2();

% This is complete (original) data records for each robot
this.LelyMap = lely_map;
this.LelyNumMap = lely_num_map;
this.Robots = robots;

this.isLelyData = true;

end