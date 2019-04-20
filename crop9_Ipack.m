function Ipack = crop9_Ipack(Ipack,wavelengthLimits)

I = Ipack.int2;

%wavelengthLimits = [rect(1),(rect(1)+rect(3))];

I2 = I(:,wavelengthLimits(1):wavelengthLimits(2));

I2 = I2-min(min(I2));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
x = 1 / max(max(I2));
I2 = I2 * x;

Ipack.crop = I2;

close all

end