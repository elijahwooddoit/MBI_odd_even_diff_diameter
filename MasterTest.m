function Ipack = MasterTest(name,wavelengthLimits)

cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating\lil folder'
status = copyfile(name, 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating');
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating'


Ipack = imageGet(name);
Ipack = crop9_Ipack(Ipack,wavelengthLimits);
Ipack = rescaleOutliers_crop(Ipack);
Ipack = smoothingFilter(Ipack);
Ipack = ridge_and_center_find(Ipack);
Ipack = find_edges(Ipack);

end