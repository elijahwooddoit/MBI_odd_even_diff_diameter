function sp = loadSPE(filename)

%
% sp=loadSPE(filename)
%   Loads PrincetonInstrument Lightfield .spe files into a matlab structure that contains
%   1D or 2D data in the sp.int(intensity) field
%   1D: 'sp.int' returns a 1D intensity array.
%   2D: 'sp.int' returns a 2D intensity matrix.
%   3D: 'sp.NumFrames' returns the frames number of a video,
%       'sp.int' has a third dim of frames, haven't tested.
%   In addition, 
%   'sp.wavelength' returns the wavelength array (column array)
%   'sp.expo_time' returns the exposure time of the experiment in [second]
%   
% This code is rewritten since the .spe file format totally changed in version 3.0
% Basically it is a mixture of binary data (as 2.0 version) and newly added XML ASCII footer.
% We need to read the offset location from the binary part, which tells the starting position of the
% XML data, and then read from there.
%  
% Known issue: Temporarily do not support multiple Region of Interest. 
% Please refer to ftp://ftp.piacton.com/Public/Manuals/Princeton%20Instruments/SPE%203.0%20File%20Format%20Specification.pdf 
% to learn how to extract the multiple data blocks.

% Author: Zhaorong Wang
% Date: 2016-02-19

%------------ Modification History --------------------
% 2015-02-22 by Zhaorong: Revised the data fields; Added the back-compatibility of SPE 2.0
% 2015-02-21 by Zhaorong: Finished the basic functionalities.

%============ Program Begins =============================================

% Test filename existence
if exist(filename)~=2 
    if exist([filename '.spe'])==2  % Append with .spe and try again
        filename=[filename '.spe'];
    else
        error(['File "' filename '" does not exist!']); 
        sp=[];
        return
    end
end

% If file exists, open it with a file id.
fid = fopen(filename, 'r');
if fid==-1 
    errordlg('Can''t read file (insufficient rights?)!');
    sp=[];
    return
end

% Read basic information from the binary header
header = readheader(fid); % readheader is a function defined later
sp = struct();
if (header.xdim == 0 || header.ydim == 0)
    errordlg('Zero data dimensions! Probably due to multiple Region of Interest chosen');
    return
end
sp.xdim = header.xdim;
sp.ydim = header.ydim;

% Read actual data and store columnwise in .int field
fseek(fid,4100,'bof');  % Move the file pointer to the beginning of the intensity data

switch header.Datatype
    case 0
        typ='float32'; 
    case 1
        typ='int32'; %seems to work for many of our files
    case 2
        typ='int16';
    case 3
        typ='uint16';
end

for kk=1:header.NumFrames
    for ii=1:header.ydim    
        sp.int(ii,:,kk)=fread(fid,header.xdim,[typ '=>double']);
    end
end

% If version>=3.0, get calibration data from the XML footer
if header.version>=3.0
    
    % Find and save the XML footer into a temp.xml file
    fseek(fid, header.XML_offset, 'bof');
    segsize = 10000;    % Read 10kb a time
    XML_fileID = fopen('temp.xml','w');
    while ~feof(fid)
        currData = fread(fid, segsize);
        if ~isempty(currData)
            fwrite(XML_fileID, currData);
        end
    end
    fclose(XML_fileID);
    
    % Search for the 'Wavelength' and 'ExposureTime' tags and save their data
    DOMnode = xmlread('temp.xml');
    sp.wavelength = readXMLdata(DOMnode, 'Wavelength');
    sp.expo_time = readXMLdata(DOMnode, 'ExposureTime')./1000; % Convert from ms to second
    % If you want to save other data important for your experiement, you can call this function
    
    if exist('temp.xml', 'file')==2 % Delete the temp file when done.
        delete('temp.xml');
    end

elseif header.version>=2.0  % For SPE2.0, read calibration data from the header
    % Create wavelength axis from calibration polynomial
    sp.wavelength=(header.Calibpoly(1)+...
          header.Calibpoly(2)*(header.startx:header.startx+header.xdim-1)+...
          header.Calibpoly(3)*(header.startx:header.startx+header.xdim-1).^2)';
    sp.expo_time = header.Exposure;
else
    errordlg('Unsupported Version! Only support 2.x and 3.x');
end

fclose(fid);

end