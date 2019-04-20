function [intersect,diameter_pixels] = diameterIntersections(ridge1,downridge,upridge,flatridge,polyDeg)

x=ridge1(:,1);
y=ridge1(:,2);

xd = downridge(:,1);
yd = downridge(:,2);
xf = flatridge(:,1);
yf = flatridge(:,2);
xu = upridge(:,1);
yu = upridge(:,2);

temp1 = num2str(polyDeg);
temp2 = 'poly';
fit_type = strcat(temp2,temp1);

[xData, yData] = prepareCurveData( xd, yd );
ft = fittype( fit_type );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
[downfit, ~] = fit( xData, yData, ft, opts );
Xd = min(xd):0.2:min(xf);
Yd = downfit(Xd);

[xData, yData] = prepareCurveData( xf, yf );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
[flatfit, ~] = fit( xData, yData, ft, opts );
Xf = min(xd):0.2:max(xu);
Yf = flatfit(Xf);
Yf_avg = ones(length(Yf),1) * mean(yf);

[xData, yData] = prepareCurveData( xu, yu );
ft = fittype( fit_type );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
[upfit, ~] = fit( xData, yData, ft, opts );
Xu = max(xf):0.2:max(xu);
Yu = upfit(Xu);

%moved to try catch loop
%[x0,y0] = intersections(Xd,Yd,Xf,Yf);
%intersect1 = [x0,y0];
%[x0,y0] = intersections(Xf,Yf,Xu,Yu);
%intersect2 = [x0,y0];

try
    [x0,y0] = intersections(Xd,Yd,Xf,Yf);
    intersect1 = [x0,y0];
    [x0,y0] = intersections(Xf,Yf,Xu,Yu);
    intersect2 = [x0,y0];
    assert(length(intersect1) > 0);
    assert(length(intersect2) > 0);
catch
    disp('Warning: could not fit fringe flat')
    [x0,y0] = intersections(Xd,Yd,Xf,Yf_avg);
    intersect1 = [x0,y0];
    [x0,y0] = intersections(Xf,Yf_avg,Xu,Yu);
    intersect2 = [x0,y0];
    assert(length(intersect1) > 0);
    assert(length(intersect2) > 0);
end

%figure; plot(flatridge(:,1),flatridge(:,2),upridge(:,1),upridge(:,2),downridge(:,1),downridge(:,2),Xd,Yd,Xf,Yf,Xu,Yu)
figure; plot(x,y,Xd,Yd,Xf,Yf,Xu,Yu,intersect1(1),intersect1(2),'o',intersect2(1),intersect2(2),'o')

intersect = [intersect1; intersect2];
diameter = intersect2(1)-intersect1(1);

%intersections1 = intersect;
diameter_pixels = diameter;

end









