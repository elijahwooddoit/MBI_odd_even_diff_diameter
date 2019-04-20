function Ipack = MasterCrop(name)

%name = input('Please enter .SPE or .TIF filename, e.g. "myFile.SPE": ','s');        % asks user for filename with file extension, must be SPE or single tiff (type uint8 has been tested)
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating\lil folder'
status = copyfile(name, 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating');
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating'


Ipack = imageGet(name);
Ipack = crop8_Ipack(Ipack);

end