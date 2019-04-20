%take a super smoothed surface
I = Ipack.crop_smoothed2;

%make a 3d plot of the image
[x,y]=size(I);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
figure;surf(xx,yy,I);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%take a slice of the image (row)
intensity_spectrum = I(391,:);

%find the first and second derivatives of the slice and appropriately
%filter them
deriv1 = gradient(intensity_spectrum);
deriv1_f = sgolayfilt(deriv1,2,21);
%figure;plot(deriv1)    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(deriv1_f)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv2 = gradient(deriv1_f);
deriv2_f = sgolayfilt(deriv2,2,21);
%figure;plot(deriv2)    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(deriv2_f)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the low areas of the second derivative plot
threshold = min(deriv2_f)*0.3;
deriv2_low = find(deriv2_f < threshold);

%find negative peaks in the low areas of the 2nd derivative plot and create
%the fringe array
deriv2_inv = deriv2_f * -1;
[~,locs] = findpeaks(deriv2_inv);
neg_pks = intersect(deriv2_low,locs);
fringe_array = zeros(length(neg_pks),10);
fringe_array(:,1) = neg_pks';
%fringe array: [alpha(1) beta(2) a(3) b(4) c(5) d(6) e(7) f(8) i(9) j(10)]
%alpha(1) is 2nd derivative peak (left peak if doublet)
%beta(2) is right 2nd derivative peak right peak if doublet or 0 if singlet
%a(3) is 1 if doublet, 0 if singlet
%b(4) left max curvature edge
%c(5) right max curvature edge
%d(6) doublet center or 0
%e(7) left zero curvature edge
%f(8) right zero curvature edge
%i(9) left ridge
%j(10) right ridge if doublet or 0 if singlet

%algorithm for finding whether a peak is a singlet or a doublet: if
%doublet, there should only be one positive peak between two adjacent
%negative peaks
[~,locs] = findpeaks(deriv2_f); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%should try this using a superfiltered version od deriv2
skip = 0;
eliminate_rows = [];
for n = 1:(length(neg_pks)-1)
    if skip == 1
        skip = 0;
        continue
    end
    span = neg_pks(n):neg_pks(n+1);
    pos_pks = intersect(locs,span);
    if length(pos_pks) == 1
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    end
end
fringe_array(eliminate_rows,:) = [];


%check for viable fringes: if both outside edges end in a peak above zero,
%then the peak is viable - eliminate unviable fringes
eliminate_rows = [];
[a,~] = size(fringe_array);
[~,locs] = findpeaks(deriv2_f,'MinPeakHeight',0.0001);
for n = 1:a
    if fringe_array(n,3) == 0
        ind = fringe_array(n,1);
        pks1 = locs(locs < ind);
        pks2 = locs(locs > ind);
        if length(pks1) < 1
            eliminate_rows = [eliminate_rows,n];
        elseif length(pks2) < 1
            eliminate_rows = [eliminate_rows,n];
        end
    elseif fringe_array(n,3) == 1
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,2);
        pks1 = locs(locs < ind1);
        pks2 = locs(locs > ind2);
        if length(pks1) < 1
            eliminate_rows = [eliminate_rows,n];
        elseif length(pks2) < 1
            eliminate_rows = [eliminate_rows,n];
        end
    end
end
fringe_array(eliminate_rows,:) = [];

%check that peak points are actually peak points by checking to make sure
%they're not on a slope and checking their height against the average
%height of the pixels they're surrounded by

%find edges of max curvature by finding the outermost positive edge peaks
%of each fringe in the 2nd derivative plot, the area between these edges
%also doubles as the effective width of the full fringe
[a,~] = size(fringe_array);
[inty,locs] = findpeaks(deriv2_f,'MinPeakHeight',0.0001);
for n = 1:a
    if fringe_array(n,3) == 0
        ind = fringe_array(n,1);
        pks1 = locs(locs < ind);
        ht1 = inty(locs < ind);
        pks2 = locs(locs > ind);
        ht2 = inty(locs > ind);
        fringe_array(n,4) = pks1(length(pks1));
        fringe_array(n,5) = pks2(1);
    elseif fringe_array(n,3) == 1
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,2);
        pks1 = locs(locs < ind1);
        ht1 = inty(locs < ind1);
        pks2 = locs(locs > ind2);
        ht2 = inty(locs > ind2);
        fringe_array(n,4) = pks1(length(pks1));
        fringe_array(n,5) = pks2(1);
    end
end

%find the center of the doublets by finding the center peak in the 2nd
%derivative fringe plot
[a,~] = size(fringe_array);
[~,locs] = findpeaks(deriv2_f);
for n = 1:a
    if fringe_array(n,3) == 1
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,2);
        range = ind1:ind2;
        peak = intersect(range,locs);
        fringe_array(n,6) = peak;
    end
end

%find the edges by max slope/zero curvature by finding the first intersect 
%of the 2nd derivative with zero (going outward from the 2nd deriv neg peak
xaxis = zeros(1,length(deriv2_f));
x = 1:length(xaxis);
%figure;plot(x,deriv2_f,x,xaxis)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x0,~,~,~] = intersections(x,deriv2_f,x,xaxis);
locs = round(x0);
[a,~] = size(fringe_array);

for n = 1:a
    if fringe_array(n,3) == 0
        ind = fringe_array(n,1);
        pks1 = locs(locs < ind);
        fringe_array(n,7) = pks1(length(pks1));
        pks2 = locs(locs > ind);
        fringe_array(n,8) = pks2(1);
    elseif fringe_array(n,3) == 1
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,2);
        pks1 = locs(locs < ind1);
        fringe_array(n,7) = pks1(length(pks1));
        pks2 = locs(locs > ind2);
        fringe_array(n,8) = pks2(1);
    end
end

%find the ridges by zero slope by taking the absolute value of the first
%derivative, flipping it, and boosting the area with the peaks
deriv2_low = find(deriv2_f < threshold);
deriv1_absval = abs(deriv1_f);
deriv1_absval = -1 * deriv1_absval;
deriv1_absval = sgolayfilt(deriv1_absval,2,5);
deriv1_absval = deriv1_absval - min(deriv1_absval);
deriv1_absval(deriv2_low) = deriv1_absval(deriv2_low) * 10;
figure;plot(deriv1_absval) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = max(deriv1_absval) * 0.5;
[~,locs] = findpeaks(deriv1_absval,'MinPeakHeight',threshold);
[a,~] = size(fringe_array);
for n = 1:a
    if fringe_array(n,3) == 0
        ind = fringe_array(n,4);
        pks = locs(locs > ind);
        fringe_array(n,9) = pks(1);
    elseif fringe_array(n,3) == 1
        ind1 = fringe_array(n,4);
        ind2 = fringe_array(n,5);
        pks1 = locs(locs > ind1);
        pks2 = locs(locs < ind2);
        fringe_array(n,9) = pks1(1);
        fringe_array(n,10) = pks2(length(pks2));
    end
end




