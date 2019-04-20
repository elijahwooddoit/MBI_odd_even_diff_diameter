function fringe_array = fringe_finder(Image,row)

%take a super smoothed surface
%I = Ipack.crop_smoothed2;
I = Image;

%make a 3d plot of the image
[x,y]=size(I);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
%figure;surf(xx,yy,I);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%take a slice of the image (row)
intensity_spectrum = I(row,:);
%figure;plot(intensity_spectrum)

%find the first and second derivatives of the slice and appropriately
%filter them
deriv1 = gradient(intensity_spectrum);
deriv1_f = sgolayfilt(deriv1,2,9);
deriv1_ff = sgolayfilt(deriv1,2,21);
%figure;plot(deriv1);title('deriv1');    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv1_f);title('deriv1 f');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv1_ff);title('deriv1 ff');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv2 = gradient(deriv1_ff);
deriv2_f = sgolayfilt(deriv2,2,9);
deriv2_ff = sgolayfilt(deriv2,2,31);
%deriv2_f(1) =.0002;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%just playing around%%%%%%
%deriv2_f(2) = .00025;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%just playing around%%%%%
%figure;plot(deriv2);title('deriv2');    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv2_f);title('deriv2 f');   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;plot(deriv2_ff);title('deriv2 ff');  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if (peak1 - 150) >0
        indexL = peak1-150;
    else
        indexL = 1;
    end
    if (peak1 + 150) < (length(intensity_spectrum) + 1)
        indexR = peak1 + 150;
    else
        indexR = length(intensity_spectrum);
    end
    temp_Im = intensity_spectrum(:,indexL:indexR);
    mean_int = mean(mean(temp_Im));
    threshold = mean_int + 0.1;
    %disp(n);disp(threshold)
    %threshold_stuff=[threshold_stuff;row,threshold];
    if int1 < threshold
        eliminate_pks = [eliminate_pks,n];
    end
end
neg_pks(eliminate_pks) = [];

fringe_array = zeros(length(neg_pks),16);   %code to establish the fringe array, also establishes fringe array 1
fringe_array(:,1) = neg_pks';
%fringe array: [alpha(1) beta(2) a(3) b(4) c(5) d(6) e(7) f(8) g(9) h(10) gamma(11) i(12) j(13) k(14) l(15) m(16)]
%alpha(1) is 2nd derivative peak (left-most peak if not singlet)
%beta(2) is 0 if singlet, second left-most peak if not singlet
%a(3) fringe classification: 0 singlet, 1 doublet, 2 for more than a doublet
%b(4) left max curvature edge
%c(5) right max curvature edge
%d(6) fringe center if not 0
%e(7) left zero curvature edge
%f(8) right zero curvature edge
%g(9) left-most ridge
%h(10) 0 if single ridge, next peak to the right if not
%gamma(11) zero if singlet, right-most neg peak index if not
%i(12) 0 if single ridge, next peak to the right if not
%j(13) 0 if single ridge, next peak to the right if not
%k(14) if more than 2 ridges, the indices of the 1st of the two highest ridges
%l(15) if more than 2 ridges, the indices of the 2nd of the two highest ridges
%m(16) number of ridges within the fringe

%Classify neg peaks and begin grouping them:
%
%classification system: 0=singlet, 1=doublet, 2=more than doublet
%
%algorithm for finding whether a peak is a singlet or a doublet: if
%doublet, there should only be one positive peak between two adjacent
%negative peaks
%
%[~,locs1] = findpeaks(deriv2_ff); eliminated part that grouped peaks using
%super filtered second derivative to avoid ordering issues
[~,locs2] = findpeaks(deriv2_f);
skip = 0;
eliminate_rows = [];
[a,~] = size(fringe_array);
for n=1:(a-1)   %make initial combine loop using neg_pks
    if skip == 1    %skip analyzing fringe array rows that are set to be combined
        skip = 0;
        continue
    end
    span = neg_pks(n):neg_pks(n+1);
    pos_pks2 = intersect(locs2,span);
    if length(pos_pks2) == 1        %code to identify (0&1) singlet vs doublet peaks using lightly smoothed 2nd deriv, also fills in fringe_array 2,3
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        %fringe_array(n,6) = pos_pks2;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    end
end
fringe_array(eliminate_rows,:) = [];
skip = 0;
combo_made = 1;
eliminate_rows = [];
while combo_made == 1   %make second combination loop using only fringe array
    [a,~] = size(fringe_array);
    for n=1:(a-1)
        if skip == 1    %skip analyzing fringe array rows that are set to be combined
            skip = 0;
            continue
        end
        if fringe_array(n,3) == 0   %determine classification of row and 2nd deriv pks and the space between them using points of max curvature
            pk_right = fringe_array(n,1);
        elseif fringe_array(n,3) == 1
            pk_right = fringe_array(n,2);
        elseif fringe_array(n,3) == 2
            pk_right = fringe_array(n,11);
        end
        pk_left = fringe_array((n+1),1);
        span = pk_right:pk_left;
        combo = 0;
        pos_pks2 = intersect(locs2,span);
        if length(pos_pks2) == 1
            combo = 1;
            skip = 1;
        end
        if combo == 1   %combine fringes that meet criteria based on their classifications
            f1 = n;
            f2 = n+1;
            class1 = fringe_array(f1,3);
            class2 = fringe_array(f2,3);
            eliminate_rows = [eliminate_rows,f2];
           if class1 == 0 && class2 == 0
                fringe_array(f1,3) = 1;
                fringe_array(f1,2) = fringe_array(f2,1);
            elseif class1 == 0 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 0 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,11);
           elseif class1 == 1 && class2 == 0
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,1);
            elseif class1 == 1 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 1 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,11);
            elseif class1 == 2 && class2 == 0
                fringe_array(f1,11) = fringe_array(f2,1);
            elseif class1 == 2 && class2 == 1
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 2 && class2 == 2
                fringe_array(f1,11) = fringe_array(f2,11);
            end
        end
    end
    fringe_array(eliminate_rows,:) = [];
    if eliminate_rows > 0
        combo_made = 1;
    else
        combo_made = 0;
    end
    eliminate_rows = [];
end

%check for viable fringes: if both outside edges end in a peak above zero,
%then the peak is viable - eliminate unviable fringes. Also fills in fringe
%array 4,5
%
%deriv2_med = sgolayfilt(deriv2,2,9);  % deriv2 f is the same thing and it's more 'standardized'
deriv2_med = deriv2_f;
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
    elseif fringe_array(n,3) == 1
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
    elseif fringe_array(n,3) == 2
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

%check instensity spectrum for places where the first derivative is zero,
%the intensity is high enough, and the high point doesn't interesct with
%any previously established fringes to try and find new fringes
xaxis = zeros(1,length(deriv1_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_f,x,xaxis);
locs = round(x0);
[a,b] = size(fringe_array);
span1 = x(deriv2_f < 0);
potential_peaks = intersect(locs,span1);
eliminate_rows = [];
for n = 1:length(potential_peaks)   %check for zero slope points with high curvature
    idx = potential_peaks(n);
    int = intensity_spectrum(idx);
    if (idx - 150) > 0
        indexL = idx-150;
    else
        indexL = 1;
    end
    if (idx + 150) < (length(intensity_spectrum) + 1)
        indexR = idx + 150;
    else
        indexR = length(intensity_spectrum);
    end
    temp_Im = intensity_spectrum(:,indexL:indexR);
    mean_int = mean(mean(temp_Im));
    threshold = mean_int + 0.1;
    %disp(n);disp(threshold)
    %threshold_stuff=[threshold_stuff;row,threshold];
    if int < threshold
        eliminate_rows = [eliminate_rows,n];
    end
end
potential_peaks(eliminate_rows,:) = [];
%
span2 = []; %eliminate the options that overlap with previously found fringes
for n = 1:a
    fringe_area = fringe_array(n,4):fringe_array(n,5);
    span2 = [span2,fringe_area];
end
bad_positions = intersect(span2,potential_peaks);
potential_peaks = setdiff(potential_peaks,bad_positions);
fringe_array2 = zeros(length(potential_peaks),b);
column_3 = zeros(length(potential_peaks),1);
[~,locs] = findpeaks(deriv2_f,'MinPeakHeight',0);
fringe_array2(:,1) = potential_peaks;
fringe_array2(:,3) = column_3;
eliminate_rows = [];
if potential_peaks > 0
    for n = 1:length(potential_peaks)
        fringe_array2(n,1) = potential_peaks(n);
        %finally add them to the main fringe array and sort them
        locsL = locs(locs < fringe_array2(n,1));
        locsR = locs(locs > fringe_array2(n,1));
        if length(locsL) > 0 && length(locsR) > 0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Operands to the || and && operators must be convertible to logical scalar values.
            fringe_array2(n,4) = locsL(length(locsL));
            fringe_array2(n,5) = locsR(1);
        else
            eliminate_rows = [eliminate_rows,n];
        end
    end
end
fringe_array2(eliminate_rows,:) = [];
%
%add fringe array 2 to original fringe array and sort based on 1st neg
%fringes
fringe_array = [fringe_array;fringe_array2];
fringe_array = sortrows(fringe_array);

%check that peak points are actually peak points by checking to make sure
%they're not on a slope
[a,~] = size(fringe_array);
eliminate_rows = 1:a;
cancel_eliminate = [];
for n = 1:a
    if fringe_array(n,3) == 0
        if intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,1)) && intensity_spectrum(fringe_array(n,5)) < intensity_spectrum(fringe_array(n,1))
            cancel_eliminate = [cancel_eliminate,n];
        end
    elseif fringe_array(n,3) == 1
        if intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,1)) && intensity_spectrum(fringe_array(n,5)) < intensity_spectrum(fringe_array(n,1))
            cancel_eliminate = [cancel_eliminate,n];
        elseif intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,2)) && intensity_spectrum(fringe_array(n,5)) < intensity_spectrum(fringe_array(n,2))
            cancel_eliminate = [cancel_eliminate,n];
        end
    elseif fringe_array(n,3) == 2
        if intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,1)) && intensity_spectrum(fringe_array(n,5)) > intensity_spectrum(fringe_array(n,1))
            cancel_eliminate = [cancel_eliminate,n];
        elseif intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,2)) && intensity_spectrum(fringe_array(n,5)) < intensity_spectrum(fringe_array(n,2))
            cancel_eliminate = [cancel_eliminate,n];
        elseif intensity_spectrum(fringe_array(n,4)) < intensity_spectrum(fringe_array(n,11)) && intensity_spectrum(fringe_array(n,5)) < intensity_spectrum(fringe_array(n,11))
        end
    end
end
eliminate_rows(cancel_eliminate) = [];
fringe_array(eliminate_rows,:) = [];

%combine peaks if the intensity between them never drops below a threshold
%
%put this in a while loop that keeps looking and eliminating until no new
%combinations are made
%
combo_made = 1;
eliminate_rows = [];
while combo_made == 1
    [a,~] = size(fringe_array);
    for n=1:(a-1)
        if skip == 1    %skip analyzing fringe array rows that are set to be combined
            skip = 0;
            continue
        end
        if fringe_array(n,3) == 0   %determine ridges-ish and the space between them using points of max curvature
            pk_right = fringe_array(n,1);
        elseif fringe_array(n,3) == 1
            pk_right = fringe_array(n,2);
        elseif fringe_array(n,3) == 2
            pk_right = fringe_array(n,11);
        end
        pk_left = fringe_array((n+1),1);
        span = pk_right:pk_left;
        inten_span = intensity_spectrum(span);
        combo = 0;
        below_threshold = [];
        m = 1;
        while 1 > length(below_threshold) && m <= length(span)  %check to see if at any point between adjacent ridges the intensity drops below a threshold
            idx = span(m);
            int = inten_span(m);
            if (idx - 150) > 0
                indexL = idx-150;
            else
                indexL = 1;
            end
            if (idx + 150) < (length(intensity_spectrum) + 1)
                indexR = idx + 150;
            else
                indexR = length(intensity_spectrum);
            end
            temp_Im = intensity_spectrum(:,indexL:indexR);
            mean_int = mean(mean(temp_Im));
            threshold = mean_int + 0.1;
            %disp(n);disp(threshold)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %threshold_stuff=[threshold_stuff;row,threshold];
            if int < threshold
                below_threshold = idx;
            end
            m = m+1;
        end
        if length(below_threshold) < 1
            combo = 1;
            skip = 1;
        end
        if combo == 1   %combine fringes that meet criteria based on their classifications
            f1 = n;
            f2 = n+1;
            class1 = fringe_array(f1,3);
            class2 = fringe_array(f2,3);
            eliminate_rows = [eliminate_rows,f2];
           if class1 == 0 && class2 == 0
                fringe_array(f1,3) = 1;
                fringe_array(f1,2) = fringe_array(f2,1);
            elseif class1 == 0 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 0 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,11);
           elseif class1 == 1 && class2 == 0
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,1);
            elseif class1 == 1 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 1 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,11);
            elseif class1 == 2 && class2 == 0
                fringe_array(f1,11) = fringe_array(f2,1);
            elseif class1 == 2 && class2 == 1
                fringe_array(f1,11) = fringe_array(f2,2);
            elseif class1 == 2 && class2 == 2
                fringe_array(f1,11) = fringe_array(f2,11);
            end
        end
    end
    fringe_array(eliminate_rows,:) = [];
    if eliminate_rows > 0
        combo_made = 1;
    else
        combo_made = 0;
    end
    eliminate_rows = [];
end

%combine all peaks that have the same edges of max curvature
%
combo_made = 1;
combo = 1;
eliminate_rows = [];
while combo_made == 1
    [a,~] = size(fringe_array);
    for n=1:(a-1)
        if skip == 1    %skip analyzing fringe array rows that are set to be combined
            skip = 0;
            continue
        end
        left1 = fringe_array(n,4);
        right1 = fringe_array(n,5);
        left2 = fringe_array((n+1),4);
        right2 = fringe_array((n+1),5);
        if left1 == left2 && right1 == right2
            combo = 1;
        else
            combo = 0;
        end
        if combo == 1   %combine fringes that meet combination criteria based on their classifications
            f1 = n;
            f2 = n+1;
            class1 = fringe_array(f1,3);
            class2 = fringe_array(f2,3);
            eliminate_rows = [eliminate_rows,f2];
            skip = 1;
            if class1 == 0 && class2 == 0
                fringe_array(row(1),3) = 1;
                fringe_array(row(1),2) = fringe_array(row(m),1);
            elseif class1 == 0 && class2 == 1
                fringe_array(row(1),3) = 2;
                fringe_array(row(1),2) = fringe_array(row(m),1);
                fringe_array(row(1),11) = fringe_array(row(m),2);
            elseif class1 == 0 && class2 == 2
                fringe_array(row(1),3) = 2;
                fringe_array(row(1),2) = fringe_array(row(m),1);
                fringe_array(row(1),11) = fringe_array(row(m),11);
            elseif class1 == 1 && class2 == 0
                fringe_array(row(1),3) = 2;
                fringe_array(row(1),11) = fringe_array(row(m),1);
            elseif class1 == 1 && class2 == 1
                fringe_array(row(1),3) = 2;
                fringe_array(row(1),11) = fringe_array(row(m),2);
            elseif class1 == 1 && class2 == 2
                fringe_array(row(1),3) = 2;
                fringe_array(row(1),11) = fringe_array(row(m),11);
            elseif class1 == 2 && class2 == 0
                fringe_array(row(1),11) = fringe_array(row(m),1);
            elseif class1 == 2 && class2 == 1
                fringe_array(row(1),11) = fringe_array(row(m),2);
            elseif class1 == 2 && class2 == 2
                fringe_array(row(1),11) = fringe_array(row(m),11);
            end
        end
    end
    fringe_array(eliminate_rows,:) = [];
    if eliminate_rows > 0
        combo_made = 1;
    else
        combo_made = 0;
    end
    eliminate_rows = [];
end

%check to make sure that each fringe has at least one zero slope point with
%negative curvature (a ridge) within its bounds 
[a,~] = size(fringe_array);
eliminate_rows = [];
xaxis = zeros(1,length(deriv1_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_f,x,xaxis);
locs = round(x0);
for n = 1:a
    boundL = fringe_array(n,4);
    boundR = fringe_array(n,5);
    span = boundL:boundR;
    possible_ridges = intersect(span,locs);
    curvatures = deriv2_f(possible_ridges);
    min_curv = min(curvatures);
    if length(possible_ridges) < 1 || min_curv > 0
        eliminate_rows = [eliminate_rows,n];
    end
end
fringe_array(eliminate_rows,:) = [];

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

%combine peaks based on the number of times the first derivative crosses
%zero between their edges of zero curvature
%
xaxis = zeros(1,length(deriv1_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_f,x,xaxis);
locs = round(x0);
combo_made = 1;
combo = 1;
eliminate_rows = [];
while combo_made == 1
    [a,~] = size(fringe_array);
    for n=1:(a-1)
        if skip == 1    %skip analyzing fringe array rows that are set to be combined
            skip = 0;
            continue
        end
        pk_right = fringe_array(n,8);
        pk_left = fringe_array((n+1),7);
        span = pk_right:pk_left;
        crossings = intersect(span,locs);
        if length(crossings) == 1  %check for single zero crossings of the first derivative
            combo = 1;
        else
            combo = 0;
        end
        if combo == 1   %check to make sure crossing is above background intensity
            idx = crossings;
            int = intensity_spectrum(crossings);
            if (idx - 150) > 0
                indexL = idx-150;
            else
                indexL = 1;
            end
            if (idx + 150) < (length(intensity_spectrum) + 1)
                indexR = idx + 150;
            else
                indexR = length(intensity_spectrum);
            end
            temp_Im = intensity_spectrum(:,indexL:indexR);
            mean_int = mean(mean(temp_Im));
            threshold = mean_int + 0.1;
            %disp(n);disp(threshold)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %threshold_stuff=[threshold_stuff;row,threshold];
            if int < threshold
                combo = 0;
            end
        end
        if combo == 1   %combine fringes that meet combination criteria based on their classifications
            f1 = n;
            f2 = n+1;
            class1 = fringe_array(f1,3);
            class2 = fringe_array(f2,3);
            eliminate_rows = [eliminate_rows,f2];
            skip = 1;
            if class1 == 0 && class2 == 0
                fringe_array(f1,3) = 1;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 0 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,2);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 0 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,2) = fringe_array(f2,1);
                fringe_array(f1,11) = fringe_array(f2,11);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 1 && class2 == 0
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,1);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 1 && class2 == 1
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,2);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 1 && class2 == 2
                fringe_array(f1,3) = 2;
                fringe_array(f1,11) = fringe_array(f2,11);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 2 && class2 == 0
                fringe_array(f1,11) = fringe_array(f2,1);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 2 && class2 == 1
                fringe_array(f1,11) = fringe_array(f2,2);
                fringe_array(f1,8) = fringe_array(f2,8);
            elseif class1 == 2 && class2 == 2
                fringe_array(f1,11) = fringe_array(f2,11);
                fringe_array(f1,8) = fringe_array(f2,8);
            end
        end
    end
    fringe_array(eliminate_rows,:) = [];
    if eliminate_rows > 0
        combo_made = 1;
    else
        combo_made = 0;
    end
    eliminate_rows = [];
end

%find ridge points by finding 1st deriv zero crossings that also have
%negative curvature, then put all that into fringe array
[a,~] = size(fringe_array);
xaxis = zeros(1,length(deriv1_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_f,x,xaxis);
locs = round(x0);
neg_curv = x(deriv2_f < 0);
eliminate_rows = [];
for n = 1:a
    pk_left = fringe_array(n,4);
    pk_right = fringe_array(n,5);
    span = pk_left:pk_right;
    possible_peaks = intersect(span,locs);
    pks = intersect(possible_peaks,neg_curv);
    pks_int = intensity_spectrum(pks);
    num_pks = length(pks);
    if num_pks > 2
        [~,idx] = maxk(pks_int,2);
        idx = sort(idx);
    end
    if num_pks == 1
        fringe_array(n,9) = pks;
        fringe_array(n,16) = 1;
    elseif num_pks == 2
        fringe_array(n,9) = pks(1);
        fringe_array(n,10) = pks(2);
        fringe_array(n,16) = 2;
    elseif num_pks == 3
        fringe_array(n,14) = idx(1);
        fringe_array(n,15) = idx(2);
        fringe_array(n,16) = 3;
        fringe_array(n,9) = pks(1);
        fringe_array(n,10) = pks(2);
        fringe_array(n,12) = pks(3);
    elseif num_pks == 4
        fringe_array(n,14) = idx(1);
        fringe_array(n,15) = idx(2);
        fringe_array(n,16) = 4;
        fringe_array(n,9) = pks(1);
        fringe_array(n,10) = pks(2);
        fringe_array(n,12) = pks(3);
        fringe_array(n,13) = pks(4);
    elseif num_pks > 4
        fringe_array(n,16) = 5;
        fringe_array(n,9) = pks(idx(1));
        fringe_array(n,10) = pks(idx(2));
    elseif num_pks == 0
        eliminate_rows = [eliminate_rows,n];
    end
end
fringe_array(eliminate_rows,:) = [];

%find the middle of the fringes with multiple peaks
[a,~] = size(fringe_array);
xaxis = zeros(1,length(deriv1_f));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_f,x,xaxis);
locs = round(x0);
pos_curv = x(deriv2_f > 0);
for n = 1:a
    pk_idx = [fringe_array(n,9),fringe_array(n,10),fringe_array(n,12),fringe_array(n,13)];
    pk_idx = pk_idx(pk_idx ~= 0);
    pk_left = pk_idx(1);
    pk_right = pk_idx(length(pk_idx));
    span = pk_left:pk_right;
    possible_peaks = intersect(span,locs);
    pks = intersect(possible_peaks,pos_curv);
    num_pks = length(pks);
    if fringe_array(n,16) == 1
        fringe_array(n,6) = fringe_array(n,9);
    elseif fringe_array(n,16) > 1
        if num_pks == 0
            fringe_array(n,6) = mean(pk_idx);
        elseif num_pks == 1
            fringe_array(n,6) = pks;
        elseif num_pks > 1
            fringe_array(n,6) = mean(pk_idx);
        end
    end
end


end