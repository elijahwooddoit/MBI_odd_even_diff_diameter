function Image_enhanced = enhanceImage(Image)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%enhance the cropped image
I1 = Image;
[y,x] = size(I1);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(X,Y);
%figure;surf(xx,yy,I1);

[xx1,yy1]=ndgrid(Y,X);
F = griddedInterpolant(xx1,yy1,I1,'cubic'); %make a cubic interpolant surface fit F
a = 1:0.1:x;
b = 1:0.1:y;
[A,B] = meshgrid(a,b);
[A1,B1] = ndgrid(b,a);
I2 = F(A1,B1);  %apply surface fit to a square grid with step sizes 1/10th as large as original
%figure;surf(A,B,I2,'linestyle','none');

I2 = I2-min(min(I2));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
o = 1 / max(max(I2));
I2 = I2 * o;

Image_rescaled = I2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%smooth the enhanced image

I = Image_rescaled;

B = imgaussfilt(I,1.5,'FilterSize',25,'FilterDomain','frequency');

I2 = smoothdata(B,1,'gaussian',50);%because each grid point is now only 1/10 original size, smoothing by 50 is equivalent to smoothing by 5 previously
I3 = smoothdata(B,2,'gaussian',50);
Ifilt1 = (I2+I3)./2;

Ifilt2 = smoothdata(Ifilt1,1,'gaussian',50);
Ifilt3 = smoothdata(Ifilt1,2,'gaussian',50);
Ifilt4 = (Ifilt2+Ifilt3)./2;

Ifilt5 = smoothdata(Ifilt4,1,'gaussian',50);
Ifilt6 = smoothdata(Ifilt4,2,'gaussian',50);
Ifilt = (Ifilt5+Ifilt6)./2;

Ifilt = Ifilt-min(min(Ifilt));
o = 1 / max(max(Ifilt));
Ifilt = Ifilt * o;

%[x,y]=size(Ifilt);
%X=1:x;
%Y=1:y;
%[xx,yy]=meshgrid(Y,X);
%figure;surf(xx,yy,Ifilt,'linestyle','none');

Image_enhanced = Ifilt;

end