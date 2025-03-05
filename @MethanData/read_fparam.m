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
p.dtype = [];
p.sheet = [];
p.amsrange = [];
p.snifrange = [];
p.head = [];
p.robid = []; % robot id

% for AMS param file
p.ams = []; % type of AMS: lely | delaval
p.device = [];
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
            case 'DTYPE'
                p.dtype = split(info{i,2},",");
                dtype = strtrim(convertCharsToStrings( p.dtype ));
                p.dtype = dtype;
            case 'SHEET'
                p.sheet = info{i,2};
            case 'AMSRANGE'
                p.amsrange = str2double( split(info{i,2},",")' );
            case 'SNIFRANGE'
                p.snifrange = str2double( split(info{i,2},",")' );
            case 'HEAD'
                p.head = split(info{i,2},",");
                head = strtrim(convertCharsToStrings( p.head ));
                p.head = head;
            case 'SNIFFER'
                p.robid = str2double( info{i,2} );
            case 'AMS'
                p.ams = string( info{i,2} );
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
            otherwise
                disp( strcat( "WARNING: not recognised parameter's KEYWORD: ", convertCharsToStrings(info{i,1}) ) );
        end
    end
end

end