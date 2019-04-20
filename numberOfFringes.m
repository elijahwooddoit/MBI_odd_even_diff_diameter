function Ipack = numberOfFringes(Ipack)

peakProminence = Ipack.peakProminence;
peakDist = Ipack.peakDist;
%
%I = Ipack.crop_smoothed;
%figure('Color','w','Name','Fringes');
%imshow(I);
%title('Select Center of Fringes');
%[~, rect] = imcrop; %crop image at fringe center, program will evaluate image using top and bottom of crops as vertical limits
%close all
%
%rect(1) = round(rect(1));
%rect(2) = round(rect(2));
%rect(3) = round(rect(3));
%rect(4) = round(rect(4));

I = Ipack.crop_smoothed2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed form smoothed to smoothed2
I2 = sum(I,1);

I2 = I2-min(min(I2));       %I2 is intensity integration of cropped file
o = 1 / max(max(I2));
I2 = I2 * o;

[pks,locs]=findpeaks(I2,'MinPeakProminence',0.2,'MinPeakDistance',10);  %find intensity peaks in I2 spectrum

x = 1:length(I2);
figure; plot(x,I2,locs,pks,'o')

fringeCenters = zeros(length(locs),1);  %find centers of fringes by looking near the intensity peaks
for n = 1:length(locs)
    center = locs(n);
    fringeCenter = fitFringeCenter(center, I, I2, Ipack);
    fringeCenters(n) = fringeCenter;
end

spectrumCenter = round(mean(fringeCenters));    %define center area as +/- 15 from avg center
lowerbound = spectrumCenter - 15;
upperbound = spectrumCenter + 15;

%I2 = I(rect(2):(rect(2)+rect(4)),:);
I3 = I(lowerbound:upperbound,:);

I3 = I3-min(min(I3));       %Ipack.int2 is the original image scaled from 0-1 intesity, base work on this file
x = 1 / max(max(I3));
I3 = I3 * x;

numFringe = zeros(1,size(I3,1)); %%%find the number of peaks within the vertical limits and interpret this as the number of fringes
numFringe = numFringe';
for w=1:length(numFringe)
    A = I3(w,:);
    [~,locs]=findpeaks(A,'MinPeakProminence',peakProminence,'MinPeakDistance',peakDist);
    if isempty(locs)
    else
        numFringe(w,:) = length(locs);
    end
end
fringeNum = mode(numFringe);

Ipack.fringeCenter = spectrumCenter;
Ipack.fringeNum = fringeNum;

end