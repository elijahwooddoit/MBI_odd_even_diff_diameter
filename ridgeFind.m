function Ipack = ridgeFind(Ipack)

peakProminence = Ipack.peakProminence;
peakDist = Ipack.peakDist;
fringeNum = Ipack.fringeNum;
Ifilt = Ipack.crop_smoothed2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed from smoothed to smoothed2

ridge = zeros(fringeNum,size(Ifilt,1)); %%%change size to number of peaks I should be detecting
ridge = ridge';
for w=1:length(ridge)
    A = Ifilt(w,:);
    [~,locs]=findpeaks(A,'MinPeakProminence',peakProminence,'MinPeakDistance',peakDist);
    if isempty(locs)
    elseif length(locs) == fringeNum
        ridge(w,:) = locs;
    end
end

[~, b] = size(ridge);
temp = 1:length(ridge); %add pixel rows to ridge matrix
temp = temp';
[c, ~] = find(ridge == 0); %finds any rows that still have dummy zero values and eliminates them
dummy = unique(c);
ridge(dummy,:) = [];
temp(dummy) = [];

a = length(temp);
ridge1 = zeros(a,3);
temp2 = zeros(a,1);
for n = 1:b
    ridge1(:,:,n) = [temp, ridge(:,n), temp2];
end

[a, ~, c] = size(ridge1);
for n = 1:a    %finds intensity values associated with all ridge points after they have been filtered down
    for m = 1:c
        temp3 = Ifilt(ridge1(n,1,m),ridge1(n,2,m));
        ridge1(n,3,m) = temp3;
    end
end

Ipack.crop_ridge = ridge1;

ridge2 = [];
for n = 1:c
    dummy2 = ridge1(:,:,n);
    ridge2 = [ridge2; dummy2];
end
 
[x,y]=size(Ifilt);  %plots ridges on 3D representation
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
figure;surf(xx,yy,Ifilt);line(ridge2(:,2),ridge2(:,1),ridge2(:,3));

end