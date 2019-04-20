%take a super smoothed surface
I = Ipack.crop_smoothed2;

%make a 3d plot of the image
[x,y]=size(I);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
figure;surf(xx,yy,I);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%take a slice of the image (row)
intensity_spectrum = I(505,:);
%figure;plot(intensity_spectrum)

%find the first and second derivatives of the slice and appropriately
%filter them
deriv1 = gradient(intensity_spectrum);
deriv1_f = sgolayfilt(deriv1,2,9);
deriv1_ff = sgolayfilt(deriv1,2,21);
%figure;plot(deriv1);title('deriv1');    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(deriv1_f);title('deriv1 f');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv1_ff);title('deriv1 ff');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv2 = gradient(deriv1_ff);
deriv2_f = sgolayfilt(deriv2,2,9);
deriv2_ff = sgolayfilt(deriv2,2,31);
deriv2_f(1) = .0002;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv2_f(2) = .00025;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv2);title('deriv2');    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(deriv2_f);title('deriv2 f');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot(deriv2_ff);title('deriv2 ff');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the low areas of the second derivative plot
threshold = min(deriv2_f)*0.1;
%threshold = -0.00075;
deriv2_low = find(deriv2_f < threshold);

%find negative peaks in the low areas of the 2nd derivative plot and create
%the fringe array, also check to make sure they're high enough to be a peak
%by taking ~base height and then adding a cutoff intensity above that
deriv2_inv = deriv2_f * -1;
[~,locs] = findpeaks(deriv2_inv);
neg_pks = intersect(deriv2_low,locs);
pks_int = intensity_spectrum(neg_pks);
eliminate_pks = [];
for n = 1:length(neg_pks)
    peak1 = neg_pks(n);
    int1 = pks_int(n);
    if (peak1 - 100) >0
        indexL = peak1-100;
    else
        indexL = 1;
    end
    if (peak1 + 100) < (length(intensity_spectrum) + 1)
        indexR = peak1 + 100;
    else
        indexR = length(intensity_spectrum);
    end
    temp_Im = I(:,indexL:indexR);
    mean_int = mean(mean(temp_Im));
    threshold = mean_int + 0.05;
    %disp(n);disp(threshold)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if int1 < threshold
        eliminate_pks = [eliminate_pks,n];
    end
end
neg_pks(eliminate_pks) = [];

fringe_array = zeros(length(neg_pks),12);
fringe_array(:,1) = neg_pks';
%fringe array: [alpha(1) beta(2) a(3) b(4) c(5) d(6) e(7) f(8) i(9) j(10) g(11) h(12)]
%alpha(1) is 2nd derivative peak (left peak if doublet)
%beta(2) is right 2nd derivative peak right peak if doublet or 0 if singlet
%a(3) fringe classification: 0 singlet, 1 doublet, 2-5 hybrid (see below)
%b(4) left max curvature edge
%c(5) right max curvature edge
%d(6) doublet center or 0
%e(7) left zero curvature edge
%f(8) right zero curvature edge
%i(9) left ridge
%j(10) right ridge if doublet or 0 if singlet
%g(11) zero if not triplet, third peak index if triplet (essentially gamma)
%h(12) which of the triplet peaks got nixed: 1,2, or 3

%algorithm for finding whether a peak is a singlet or a doublet: if
%doublet, there should only be one positive peak between two adjacent
%negative peaks
%classification system: 0=singlet, 1=doublet, 2=skew left, 3=skew right,
%4=hybrid ridge left, 5=hybrid ridge right, 6=triplet
%
%
%Classify cases 0, 1, and 6 and fill in doublet center
[~,locs1] = findpeaks(deriv2_ff);
[~,locs2] = findpeaks(deriv2_f);
skip = 0;
eliminate_rows = [];
for n = 1:(length(neg_pks)-1)
    if skip == 1
        skip = 0;
        continue
    elseif skip == 2
        skip = 1;
        continue
    end
    span = neg_pks(n):neg_pks(n+1);
    pos_pks1 = intersect(locs1,span);
    pos_pks2 = intersect(locs2,span);
    if length(pos_pks1) == 1    %code to identify (0&1) singlet vs doublet peaks, also fills in 1,2,3,6 in fringe array
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        fringe_array(n,6) = pos_pks1;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    elseif length(pos_pks2) == 1
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        fringe_array(n,6) = pos_pks2;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    end
    if n+2 <= length(neg_pks) && fringe_array(n,3) == 1   %code to indentify (6) triplet peaks, also fills in 1,2,3,6,9,10,11,12 in fringe array
        span1 = neg_pks(n+1):neg_pks(n+2);
        pos_pks1 = intersect(locs1,span1);
        pos_pks2 = intersect(locs2,span1);
        if length(pos_pks1) == 1
            fringe_array(n,3) = 6;
            pks_ind = [neg_pks(n),neg_pks(n+1),neg_pks(n+2)];
            pks_ht = intensity_spectrum(pks_ind);
            low_ind = find(pks_ht == min(pks_ht));
            elim_pk = pks_ind(low_ind);
            fringe_array(n,11) = pks_ind(3);
            fringe_array(n,12) = low_ind;
            pks = pks_ind;
            pks(low_ind) = [];
            fringe_array(n,6) = round(mean(pks(1),pks(2)));
            fringe_array(n,1) = pks_ind(1);
            fringe_array(n,2) = pks_ind(2);
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
            skip = 2;
            eliminate_rows = [eliminate_rows,(n+1),(n+2)];
            error('!!!triplet detected!!!')
        elseif length(pos_pks2) == 1
            fringe_array(n,3) = 6;
            pks_ind = [neg_pks(n),neg_pks(n+1),neg_pks(n+2)];
            pks_ht = intensity_spectrum(pks_ind);
            low_ind = find(pks_ht == min(pks_ht));
            elim_pk = pks_ind(low_ind); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pks to pks_ind
            fringe_array(n,11) = pks_ind(3);
            fringe_array(n,12) = low_ind;
            pks = pks_ind;
            pks(low_ind) = [];
            fringe_array(n,6) = round(mean(pks(1),pks(2)));
            fringe_array(n,1) = pks_ind(1);
            fringe_array(n,2) = pks_ind(2);
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
            skip = 2;
            eliminate_rows = [eliminate_rows,(n+1),(n+2)];
            %error('!!!triplet detected!!!')
        end
    end
end
fringe_array(eliminate_rows,:) = [];
%
%classify cases (2&3) skewed peaks, leave second peak empty for these cases, fills in 1,2,3,6 and one ridge of fringe array
[a,~] = size(fringe_array);
for n = 1:a
    if fringe_array(n,3) == 1
        if deriv2_f(fringe_array(n,6)) < 0
            fringe_array(n,3) = 2;
        end
    end
end
for n = 1:a
    if fringe_array(n,3) == 2
        if intensity_spectrum(fringe_array(n,2)) > intensity_spectrum(fringe_array(n,1))
            fringe_array(n,3) = 3;
            fringe_array(n,9) = 0;
        else
            fringe_array(n,10) = 0;
        end
    end
end
%
%classify cases (4&5) hybrid peaks, define minor peak as the corresponding neg curvature peak, fills in 1,2,3,6 and minor ridge of fringe array
[~,locs1] = findpeaks(deriv1_ff);
xaxis = zeros(1,length(deriv1_ff));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_ff,x,xaxis);
locs = round(x0);
for n = 1:a
    if fringe_array(n,3) == 1
        span = fringe_array(n,1):fringe_array(n,2);
        pks1 = intersect(locs1,span);
        zeros0 = intersect(locs,span);
        if length(zeros0) == 1
            fringe_array(n,3) = 4;
            if zeros0 > pks1
                fringe_array(n,3) = 5;
                fringe_array(n,9) = fringe_array(n,1);
            else
                fringe_array(n,10) = fringe_array(n,2);
            end
        end
    end
end

%check for viable fringes: if both outside edges end in a peak above zero,
%then the peak is viable - eliminate unviable fringes
%deriv2_med = sgolayfilt(deriv2,2,9);  - deriv2 f is the same thing and it's
%more 'standardized'
deriv2_med = deriv2_f;
%figure;plot(deriv2_med);title('deriv2 med');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eliminate_rows = [];
[a,~] = size(fringe_array);
[~,locs] = findpeaks(deriv2_med,'MinPeakHeight',0);
for n = 1:a
    if fringe_array(n,3) == 0
        ind = fringe_array(n,1);
        pks1 = locs(locs < ind);
        pks2 = locs(locs > ind);
        if length(pks1) < 1
            eliminate_rows = [eliminate_rows,n];
        elseif length(pks2) < 1
            eliminate_rows = [eliminate_rows,n];
        else
            fringe_array(n,4) = pks1(length(pks1));
            fringe_array(n,5) = pks2(1);
        end
    elseif fringe_array(n,3) > 0 && fringe_array(n,3) < 6
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,2);
        pks1 = locs(locs < ind1);
        pks2 = locs(locs > ind2);
        if length(pks1) < 1
            eliminate_rows = [eliminate_rows,n];
        elseif length(pks2) < 1
            eliminate_rows = [eliminate_rows,n];
        else
            fringe_array(n,4) = pks1(length(pks1));
            fringe_array(n,5) = pks2(1);
        end
    elseif fringe_array(n,3) == 6
        ind1 = fringe_array(n,1);
        ind2 = fringe_array(n,11);
        pks1 = locs(locs < ind1);
        pks2 = locs(locs > ind2);
        if length(pks1) < 1
            eliminate_rows = [eliminate_rows,n];
        elseif length(pks2) < 1
            eliminate_rows = [eliminate_rows,n];
        else
            fringe_array(n,4) = pks1(length(pks1));
            fringe_array(n,5) = pks2(1);
        end
    end
end
fringe_array(eliminate_rows,:) = [];

%find edges of max curvature by finding the outermost positive edge peaks
%of each fringe in the 2nd derivative plot, the area between these edges
%also doubles as the effective width of the full fringe - added to fringe
%viability test
%[a,~] = size(fringe_array);
%[inty,locs] = findpeaks(deriv2_f,'MinPeakHeight',0.0001);
%for n = 1:a
%    if fringe_array(n,3) == 0
%        ind = fringe_array(n,1);
%        pks1 = locs(locs < ind);
%        ht1 = inty(locs < ind);
%        pks2 = locs(locs > ind);
%        ht2 = inty(locs > ind);
%        fringe_array(n,4) = pks1(length(pks1));
%        fringe_array(n,5) = pks2(1);
%    elseif fringe_array(n,3) == 1
%        ind1 = fringe_array(n,1);
%        ind2 = fringe_array(n,2);
%        pks1 = locs(locs < ind1);
%        ht1 = inty(locs < ind1);
%        pks2 = locs(locs > ind2);
%        ht2 = inty(locs > ind2);
%        fringe_array(n,4) = pks1(length(pks1));
%        fringe_array(n,5) = pks2(1);
%    end
%end

%check that peak points are actually peak points by checking to make sure
%they're not on a slope
%eliminate_rows = [];
[a,~] = size(fringe_array);
eliminate_rows = [];
for n = 1:a
    if fringe_array(n,3) == 0
        if intensity_spectrum(fringe_array(n,4)) > intensity_spectrum(fringe_array(n,1)) || intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(fringe_array(n,1))
            eliminate_rows = [eliminate_rows,n];
        end
    elseif fringe_array(n,3) > 0 && fringe_array(n,3) < 6
        if intensity_spectrum(fringe_array(n,4)) > intensity_spectrum(fringe_array(n,1)) || intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(fringe_array(n,1))
            eliminate_rows = [eliminate_rows,n];
        elseif intensity_spectrum(fringe_array(n,4)) > intensity_spectrum(fringe_array(n,2)) || intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(fringe_array(n,2))
            eliminate_rows = [eliminate_rows,n];
        end
    elseif fringe_array(n,3) == 6
        ind_all = [fringe_array(n,1),fringe_array(n,2),fringe_array(n,3)];
        ind_all(fringe_array(n,12)) = [];
        if intensity_spectrum(fringe_array(n,4)) > intensity_spectrum(ind_all(1)) || intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(ind_all(1))
            eliminate_rows = [eliminate_rows,n];
        elseif intensity_spectrum(fringe_array(n,4)) > intensity_spectrum(ind_all(2)) || intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(ind_all(2))
            eliminate_rows = [eliminate_rows,n];
        end
    end
end
fringe_array(eliminate_rows,:) = [];

%find the center of the doublets by finding the center peak in the 2nd
%derivative fringe plot - already done when categorizing peaks
%[a,~] = size(fringe_array);
%[~,locs] = findpeaks(deriv2_f);
%for n = 1:a
%    if fringe_array(n,3) == 1
%        ind1 = fringe_array(n,1);
%        ind2 = fringe_array(n,2);
%        range = ind1:ind2;
%        peak = intersect(range,locs);
%        fringe_array(n,6) = peak;
%    end
%end

%find the edges by max slope/zero curvature by finding the first intersect 
%of the 2nd derivative with zero (going inward from the edges of max curvature)
xaxis = zeros(1,length(deriv2_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv2_f,x,xaxis);
locs = round(x0);
[a,~] = size(fringe_array);
for n = 1:a
    ind1 = fringe_array(n,4);
    ind2 = fringe_array(n,5);
    pks1 = locs(locs > ind1);
    fringe_array(n,7) = pks1(1);
    pks2 = locs(locs < ind2);
    fringe_array(n,8) = pks2(length(pks2));
end

%find the edges by max slope/zero curvature by finding the first intersect 
%of the 2nd derivative with zero (going outward from the 2nd deriv neg peak
%- revised to go inward from outer edges and be applicable to all cases 1-6
%xaxis = zeros(1,length(deriv2_f));
%x = 1:length(xaxis);
%figure;plot(x,deriv2_f,x,xaxis)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[x0,~,~,~] = intersections(x,deriv2_f,x,xaxis);
%locs = round(x0);
%[a,~] = size(fringe_array);
%for n = 1:a
%    if fringe_array(n,3) == 0
%        ind = fringe_array(n,1);
%        pks1 = locs(locs < ind);
%        fringe_array(n,7) = pks1(length(pks1));
%        pks2 = locs(locs > ind);
%        fringe_array(n,8) = pks2(1);
%    elseif fringe_array(n,3) == 1
%        ind1 = fringe_array(n,1);
%        ind2 = fringe_array(n,2);
%        pks1 = locs(locs < ind1);
%        fringe_array(n,7) = pks1(length(pks1));
%        pks2 = locs(locs > ind2);
%        fringe_array(n,8) = pks2(1);
%    end
%end

%find the ridges by zero slope by taking the absolute value of the first
%derivative, flipping it, and boosting the area with the peaks
deriv1_absval = abs(deriv1_f);
deriv1_absval = -1 * deriv1_absval;
deriv1_absval = sgolayfilt(deriv1_absval,2,5);
deriv1_absval = deriv1_absval - min(deriv1_absval);
%figure;plot(deriv1_absval) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv1_absval(deriv2_low) = deriv1_absval(deriv2_low) * 10;
figure;plot(deriv1_absval) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = max(deriv1_absval) * 0.5;
[~,locs] = findpeaks(deriv1_absval,'MinPeakHeight',threshold);
%xaxis = zeros(1,length(deriv2_f));
%x = 1:length(xaxis);
%[x0,~,~,~] = intersections(x,deriv2_f,x,xaxis);
%locs2 = round(x0);
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
    elseif fringe_array(n,3) == 2 || fringe_array(n,3) == 4
        ind = fringe_array(n,4);
        pks = locs(locs > ind);
        fringe_array(n,9) = pks(1);
    elseif fringe_array(n,3) == 3 || fringe_array(n,3) == 5
        ind = fringe_array(n,5);
        pks = locs(locs < ind);
        fringe_array(n,10) = pks(length(pks));
    elseif fringe_array(n,3) == 6
        ind1 = fringe_array(n,4);
        ind2 = fringe_array(n,5);
        span1 = ind1:ind2;
        %ridge_pts = intersect(locs2,span1);
        pks = intersect(span1,locs);
        if length(pks) == 3
            pks(fringe_array(n,12)) = [];
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
        elseif length(pks) == 2
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
        elseif length(pks) == 1
            fringe_array(n,9) = pks(1);
            fringe_array(n,3) = 1;
        end
    end
end




