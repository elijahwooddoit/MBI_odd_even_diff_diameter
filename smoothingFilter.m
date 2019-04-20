function Ipack = smoothingFilter(Ipack)

I = Ipack.crop_rescaled;
[x,y]=size(I);
X=1:x;
Y=1:y;


B = imgaussfilt(I,1.5,'FilterSize',101,'FilterDomain','frequency');

%for Ifilt
I2 = smoothdata(B,1,'gaussian',20);
I3 = smoothdata(B,2,'gaussian',20);
Ifilt1 = (I2+I3)./2;

Ifilt2 = smoothdata(Ifilt1,1,'gaussian',20);
Ifilt3 = smoothdata(Ifilt1,2,'gaussian',20);
Ifilt4 = (Ifilt2+Ifilt3)./2;

Ifilt5 = smoothdata(Ifilt4,1,'gaussian',20);
Ifilt6 = smoothdata(Ifilt4,2,'gaussian',20);
Ifilt = (Ifilt5+Ifilt6)./2;

Ifilt = imgaussfilt(Ifilt,1.5,'FilterSize',101,'FilterDomain','frequency');

%Ifilt8 = smoothdata(Ifilt7,1,'gaussian',5);
%Ifilt9 = smoothdata(Ifilt7,2,'gaussian',5);
%Ifilt10 = (Ifilt8+Ifilt9)./2;

%for Ifilt2
Ifilt8 = smoothdata(B,1,'gaussian',30);
Ifilt9 = smoothdata(B,2,'gaussian',20);
Ifilt10 = (Ifilt8+Ifilt9)./2;

Ifilt11 = smoothdata(Ifilt10,1,'gaussian',30);
Ifilt12 = smoothdata(Ifilt10,2,'gaussian',10);
Ifilt13 = (Ifilt11+Ifilt12)./2;

Ifilt14 = smoothdata(Ifilt13,1,'gaussian',30);
Ifilt15 = smoothdata(Ifilt13,2,'gaussian',10);
Ifilt2 = (Ifilt14+Ifilt15)./2;

Ifilt2 = smoothdata(Ifilt2,2,'gaussian',30);

%Ifilt = smoothdata(Ifilt10,2,'gaussian',20);

[xx,yy]=meshgrid(Y,X);
%figure;surf(xx,yy,I);
%figure;surf(xx,yy,B);
%figure;surf(xx,yy,Ifilt);

Ifilt = Ifilt-min(min(Ifilt));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
o = 1 / max(max(Ifilt));
Ifilt = Ifilt * o;

Ifilt2 = Ifilt2-min(min(Ifilt2));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
o = 1 / max(max(Ifilt2));
Ifilt2 = Ifilt2 * o;

Ipack.crop_smoothed = Ifilt;
Ipack.crop_smoothed2 = Ifilt2;

end