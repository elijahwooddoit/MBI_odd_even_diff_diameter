function Ipack = smoothingFilter2(Ipack)

I = Ipack.crop_rescaled;
[x,y]=size(I);
X=1:x;
Y=1:y;


B = imgaussfilt(I,1.5,'FilterSize',25,'FilterDomain','frequency');

I2 = smoothdata(B,1,'gaussian',5);
I3 = smoothdata(B,2,'gaussian',5);
Ifilt1 = (I2+I3)./2;

Ifilt2 = smoothdata(Ifilt1,1,'gaussian',5);
Ifilt3 = smoothdata(Ifilt1,2,'gaussian',5);
Ifilt4 = (Ifilt2+Ifilt3)./2;

Ifilt5 = smoothdata(Ifilt4,1,'gaussian',5);
Ifilt6 = smoothdata(Ifilt4,2,'gaussian',5);
Ifilt = (Ifilt5+Ifilt6)./2;

%Ifilt8 = smoothdata(Ifilt7,1,'gaussian',5);
%Ifilt9 = smoothdata(Ifilt7,2,'gaussian',5);
%Ifilt10 = (Ifilt8+Ifilt9)./2;

%Ifilt = smoothdata(Ifilt10,2,'gaussian',20);

%[xx,yy]=meshgrid(Y,X);
%figure;surf(xx,yy,B);
%figure;surf(xx,yy,Ifilt);

%Ifilt = B;

Ifilt = Ifilt-min(min(Ifilt));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
o = 1 / max(max(Ifilt));
Ifilt = Ifilt * o;

Ipack.crop_smoothed2 = Ifilt;

end