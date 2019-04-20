function Ipack = imageGet(filename)

% Reads Tiff and SPE files from same folder it is located in and stores
% them in a structure called Ipack that contains dimensional, intensity,
% wavelength, and exposure time info if available. Output intensity files
% are greyscale

%name = input('Please enter .SPE or .TIF filename, e.g. "myFile.SPE": ','s');        % asks user for filename with file extension, must be SPE or single tiff (type uint8 has been tested)
name = filename;

if length(name) == 1
    return
end

if strcmp(name(end-3:end),'.TIF') == 1 || strcmp(name(end-3:end),'.tif') == 1
    I = imread(name);       % loads the Tiff image into memory
    I2 = double(I);      % converts uint8 type variables to double type
    y = 1 / max(max(I2));
    I2 = I2 * y;     % normalizes intensity array to a 0-1 scale
    Ipack = struct;     % creates structure called Ipack and stores all relevant information in structure fields
    [Ipack.ydim, Ipack.xdim] = size(I2);
    Ipack.int = I2;
    Ipack.wavelength = 'unknown';
    Ipack.expo_time = 'unknown';
    Ipack.file_type = '.tif';
    
elseif strcmp(name(end-3:end),'.SPE') == 1 || strcmp(name(end-3:end),'.spe') == 1
    Ipack = loadSPE(name);     % loads spe file into a structure called Ipack using function by Zhaorong Wang and then converts intensity array to a 0-1 scale
    y = 1 / max(max(Ipack.int));
    Ipack.int = Ipack.int * y;
    Ipack.file_type = '.spe';

else
    disp('!!! File name incorrectly entered!!!');     % error message if one of the two conditions above aren't triggered
    Ipack = [];
    return
end
Ipack.fileName = name;

int2 = Ipack.int;
int2 = int2-min(min(int2));
x = 1 / max(max(int2));
int2 = int2 * x;
Ipack.int2=int2;        %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file

end