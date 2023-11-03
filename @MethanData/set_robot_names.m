function set_robot_names( this, names )
% Set expected robots' names.
% We need this method only in case if the expected robot names
% are totally different from the default ones;
% for example, there are names like 104, 105, ... but there are
% no names like 101, 102, and 103 present in a data.
this.expected_robots = names;
end
