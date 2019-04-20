function Ipack = crop6_Ipack(Ipack)

I = Ipack.int2;
J=imadjust(I);

figure('Color','w','Name','Fringes');       % creates an empty figure with white background and figure name
imshow(J)       % shows the loaded tiff image in the previously empty figure
title('Select Region of Interest in Image');        %adds title to figure

[~, rect] = imcrop;            
% asks the user to crop with a rectangle, where rect is defined by the
% vector [x0 y0 width height], to finish cropping double click on
% rectangle, outputs cropped image in format of each pixel occupies a
% corresponding matrix site with a single intensity value

rect(1) = round(rect(1)); %rounds all files to be accurate taking pixels as opposed to continuous numbers into account
rect(2) = round(rect(2));
rect(3) = round(rect(3));
rect(4) = round(rect(4));

I2 = I(:,rect(1):(rect(1)+rect(3)));

I2 = I2-min(min(I2));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
x = 1 / max(max(I2));
I2 = I2 * x;

Ipack.crop = I2;

close all

end