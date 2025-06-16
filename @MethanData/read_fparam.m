function p = read_fparam( this, filename )

fid = fopen( filename );

tline = '0';
msg = '0';
i = 1;
info = {};
pos = 0;

while ( tline ~= -1 )
    tline = fgetl(fid);

    while ( isempty(tline) )
        tline = fgetl(fid);
    end

    if ( ~isempty(tline) && tline(1,1) == '$' )

        pos = ftell(fid);
        msg = fgetl(fid);

        while ( isempty(msg) || msg(1,1) == '#' )
            pos = ftell(fid);
            msg = fgetl(fid);
        end
        if ( msg(1,1) == '$' )
            fseek(fid,pos,'bof');
        else
            if (msg(1,1) ~= -1)
                info{i,2} = msg;
            end
        end
        info{i,1} = tline(2:end);
    end
    i = i + 1;
end

p.amsfile = [];
p.sniffile = [];
p.amsrange = [];
p.snifrange = [];
p.robid = []; % robot id, AMS station as DevAddress
p.snif_id = []; % sniffer id, optional
p.farm_id = []; % farm id, optional

% for AMS param file
p.device = []; % col name for robot id
p.id = [];
p.time = [];

% for GEDA
p.deltaT = [];
p.siglength = [];

p.wndlength = [];
p.eigenvalues = [];
p.denoise = [];

% where to write the resulting files
p.outpath = [];

% for studying of downsampling
p.dwnsampl = [];

% type of trait
p.trait = [];

% for outliers filtering in write_traits.m
p.filtering = [];

for i = 1:size(info,1)
    if ~isempty(info{i,1}) && ~isempty(info{i,2})
        switch info{i,1}
            case 'AMSFILE'
                p.amsfile = split(info{i,2},",");%info{i,2};
                pfile = strtrim(convertCharsToStrings( p.amsfile ));
                p.amsfile = pfile;
            case 'SNIFFILE'
                p.sniffile = split(info{i,2},",");%info{i,2};
                pfile = strtrim(convertCharsToStrings( p.sniffile ));
                p.sniffile = pfile;
            case 'AMSRANGE'
                p.amsrange = str2double( split(info{i,2},",")' );
            case 'SNIFRANGE'
                p.snifrange = str2double( split(info{i,2},",")' );
            case 'SNIFID'
                p.snif_id = string( info{i,2} );
            case 'FARMID'
                p.farm_id = string( info{i,2} );
           case 'ROBOTID'
                p.robid = str2double( info{i,2} );
            case 'DEVICECOL'
                p.device = string( info{i,2} );
            case 'IDCOL'
                p.id = string( info{i,2} );
            case 'TIMECOL'
                p.time = split(info{i,2},",");
                ptime = strtrim(convertCharsToStrings( p.time ));
                p.time = ptime;
            case 'SAMPLING'
                p.deltaT = str2double( info{i,2} );
            case 'SIGLENGTH'
                p.siglength = str2double( info{i,2} );
            case 'DENOISE'
                p.denoise = str2double( info{i,2} );
            case 'SSA1'
                p.wndlength = str2double( info{i,2} );
            case 'SSA2'
                p.eigenvalues = str2double( info{i,2} );
            case 'OUTPATH'
                p.outpath = split(info{i,2},",");
                pfile = strtrim(convertCharsToStrings( p.outpath ));
                p.outpath = pfile;
            case 'DWNSAMPLING'
                p.dwnsampl = str2double( info{i,2} );
            case 'TRAIT'
                p.trait = str2double( info{i,2} );
            case 'FILTERING'
                filtering = split(info{i,2},",");
                filtering = strtrim(convertCharsToStrings( filtering ));
                p.filtering = str2double(filtering);
            otherwise
                disp( strcat( "WARNING: not recognised parameter's KEYWORD: ", convertCharsToStrings(info{i,1}) ) );
        end
    end
end

end