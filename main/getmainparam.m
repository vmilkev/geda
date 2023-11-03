function [ams_f,snif_f,deltaT,siglength] = getmainparam( filename )

fid = fopen( filename );

tline = '0';
msg = '0';
i = 1;
info = {};
pos = 0;

ams_f = [];
snif_f = [];
deltaT = [];
siglength = [];

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

for i = 1:size(info,1)
    if ( ~isempty(info{i,1}) )
        switch info{i,1}
            case 'AMSFILE'
                ams_f = info{i,2};
                pfile = strtrim(convertCharsToStrings( ams_f ));
                ams_f = pfile;
            case 'SNIFFILE'
                snif_f = info{i,2};
                pfile = strtrim(convertCharsToStrings( snif_f ));
                snif_f = pfile;
            case 'SAMPLING'
                deltaT = str2double( info{i,2} );
            case 'SIGLENGTH'
                siglength = str2double( info{i,2} );
        end
    end
end

end