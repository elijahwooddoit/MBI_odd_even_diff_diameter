function Ipack = MasterFile(name,wavelengthLimits)

%name = input('Please enter .SPE or .TIF filename, e.g. "myFile.SPE": ','s');        % asks user for filename with file extension, must be SPE or single tiff (type uint8 has been tested)
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating\lil folder'
status = copyfile(name, 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating');
cd 'C:\Users\Mike\Documents\MATLAB\Radius Measurement\Automation attempts\2019.01.11 revision for fully automating'


Ipack = imageGet(name);
Ipack = crop9_Ipack(Ipack,wavelengthLimits);
Ipack = crop7_Ipack(Ipack);
Ipack = rescaleOutliers_crop(Ipack);
Ipack = smoothingFilter(Ipack);

peakProminence = 0.1;
peakDist = 10;
Ipack.peakProminence=peakProminence;
Ipack.peakDist=peakDist;

Ipack = numberOfFringes(Ipack);
Ipack = ridgeFind(Ipack);

slopeFilter = 0.5;  %enter cutoff slope to define fringe edge
flatFilter = 0.25;  %enter cutoff slope to define flat of fringe
polynomialDegree = 1;   %enter polynomial degree for diameter estimate based on fitting / intersections
Ipack.slopeFilter=slopeFilter;
Ipack.flatFilter=flatFilter;
Ipack.polynomialDegree=polynomialDegree;

Ipack = ridges_diameter(Ipack);

Ipack = isolate3Dfringe(Ipack);
Ipack = smoothingFilter2(Ipack);
Ipack = imageEnhancer(Ipack);
Ipack = isolationRecalculation(Ipack);
Ipack = preparethefringes2(Ipack);

%Ipack = preparethefringes(Ipack);
%fringe_centerlines_multiple_molds2;
Ipack = fringe_centerlines_multiple_molds2(Ipack);
Ipack = errorPlotAnalyzer(Ipack);

delete(name)

end