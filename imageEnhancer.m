function Ipack = imageEnhancer(Ipack)

I1 = Ipack.crop_smoothed2;
[y,x] = size(I1);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(X,Y);
%figure;surf(xx,yy,I1);

[xx1,yy1]=ndgrid(Y,X);
F = griddedInterpolant(xx1,yy1,I1,'cubic');
a = 1:0.1:x;
b = 1:0.1:y;
[A,B] = meshgrid(a,b);
[A1,B1] = ndgrid(b,a);
I2 = F(A1,B1);
%figure;surf(A,B,I2);

I2 = I2-min(min(I2));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
o = 1 / max(max(I2));
I2 = I2 * o;

Ipack.ImEnhance = I2;

%%%%%%%%%%%%%%%%%%%%%%%%%% switched to griddedInterpolant Matlab function

%I1 = Ipack.crop_smoothed2;

%[a,b] = size(I1);

%c = (a-1)*10+1;
%d = (b-1)*10+1;
%I2 = zeros(a,d);
%for n = 1:a
%    y = I1(n,:);
%    x = 1:length(y);
%    [xData, yData] = prepareCurveData( x, y );
%    ft = 'pchipinterp';
%    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
%    xfit = 1:0.1:b;
%    yfit = fitresult(xfit);
%    yfit = yfit';
%    I2(n,:) = yfit;
%end
%
%I3 = zeros(c,d);
%for n = 1:d
%    y = I2(:,n);
%    y = y';
%    x = 1:length(y);
%    [xData, yData] = prepareCurveData( x, y );
%    ft = 'pchipinterp';
%    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
%    xfit = 1:0.1:a;
%    yfit = fitresult(xfit);
%    yfit = yfit';
%    I3(:,n) = yfit;
%end
%
%Ipack.ImEnhance = I3;

end



