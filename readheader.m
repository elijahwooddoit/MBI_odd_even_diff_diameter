function header = readheader(fid)

%each entry in headerinfo corresponds to data that is read into the header
%structure of the resulting Matlab spectrum structure
%fields are: Name, Offset (Byte number in .spe file), Type (Datatype),
%Length (for Arrays of Type), and Load (whether or not this item should be
%read)
c=1;
headerinfo(c)=struct('Name','Exposure', 'Offset',10,  'Type','float',  'MType','double','Length',1,'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','Date',     'Offset',20,  'Type','char',   'MType','char',  'Length',10,'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','Time',     'Offset',172, 'Type','char',   'MType','char',  'Length',7,'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','xdim',     'Offset',42,  'Type','uint16', 'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','ydim',     'Offset',656, 'Type','uint16', 'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','Datatype', 'Offset',108, 'Type','int16',  'MType','double', 'Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','NumFrames','Offset',1446,'Type','int32',  'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','XML_offset','Offset',678,'Type','uint64',  'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','version',  'Offset',1992,'Type','float32',  'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','startx',   'Offset',1512,'Type','uint16', 'MType','double','Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','Calibpoly','Offset',3263,'Type','float64','MType','double','Length',6, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','BckGrdApplied','Offset',150,'Type','char','MType','double', 'Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','FltFldApplied','Offset',706,'Type','char','MType','double', 'Length',1, 'Transpose',1,'Load',1); c=c+1;
headerinfo(c)=struct('Name','Gain',     'Offset',198,  'Type','int16', 'MType','double','Length',1, 'Transpose',1,'Load',1);

for i=1:length(headerinfo)
    if headerinfo(i).Load
        fseek(fid,headerinfo(i).Offset,'bof');
        if headerinfo(i).Transpose
            header.(headerinfo(i).Name)=fread(fid,headerinfo(i).Length,[headerinfo(i).Type '=>' headerinfo(i).MType])';
        else
            header.(headerinfo(i).Name)=fread(fid,headerinfo(i).Length,[headerinfo(i).Type '=>' headerinfo(i).MType]);
        end 
    end
end

return