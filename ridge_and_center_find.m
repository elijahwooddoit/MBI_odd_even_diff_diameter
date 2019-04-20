function Ipack = ridge_and_center_find(Ipack)

%open up the image and define its size
I = Ipack.crop_smoothed2;
[a,b] = size(I);
siz = [a b];
%
fringes = [];
for n = 1:a
    %disp(n)
    %store fringe array in fringes, which is just the collection of fringe
    %array matrices from every row
    fringe_array = fringe_finder(I,n);
    [c,~] = size(fringe_array);
    fringe_array = [fringe_array, (ones(c,1)*n)];
    fringes = [fringes; fringe_array];
end
linearInd = sub2ind(siz, fringes(:,17), fringes(:,6));
fringes = [fringes,linearInd];
%fringe array: [alpha(1) beta(2) a(3) b(4) c(5) d(6) e(7) f(8) g(9) h(10) gamma(11) i(12) j(13) k(14) l(15) m(16) n(17) o(18)]
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
%n(17) row fringe is situated in
%o(18) linear index of fringe middle

%make new blank image
[a,b] = size(I);
Ibounds = zeros(a,b);

%put all pixels within outer bounds as 100% brightness (value=1)
[c,~] = size(fringes);
for n = 1:c
    Ibounds(fringes(n,17),fringes(n,4):fringes(n,5)) = 1;
end

%filter the bounds image
B = imgaussfilt(Ibounds,1.5,'FilterSize',25,'FilterDomain','frequency');
I2 = smoothdata(B,1,'gaussian',100);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%filter size for fringe fattening
I3 = smoothdata(B,2,'gaussian',100);
Ifilt = (I2+I3)./2;

%threshold the filtered image
[a,b] = size(I);
Ibounds2 = zeros(a,b);
Ibounds2(Ifilt > 0.2) = 1;

%sort fringe rows into columns by finding connected pixels in a binary 
%image and then finding intersecting linear indices from the fringes
CC = bwconncomp(Ibounds2,4);
fringe_container = CC.PixelIdxList;
e = length(fringe_container);
fringes_sorted = cell(1,e);
[a,b] = size(I);
siz = [a,b];
for n = 1:e
    [~,ia,~] = intersect(fringes(:,18),fringe_container{n});
    fringes_sorted{n} = fringes(ia,:);
    fringes_sorted{n} = sortrows(fringes_sorted{n},17);
end

%create an intensity spectrum using the filtered image and filter the
%signal
int_spect = sum(I);
int_spect = sgolayfilt(int_spect,2,9);
int_spect = int_spect - min(int_spect);
int_spect = (1/max(int_spect)).*int_spect;

%find peaks in the intensity signal and make sure they're tall enough,
%these should correspond to the number of fringes in the image
[pks,locs] = findpeaks(int_spect);
eliminate_pks = [];
for n = 1:length(locs)
    peak1 = locs(n);
    int1 = pks(n);
    if (peak1 - 100) > 0
        indexL = peak1-100;
    else
        indexL = 1;
    end
    if (peak1 + 100) <= length(int_spect)
        indexR = peak1 + 100;
    else
        indexR = length(int_spect);
    end
    temp_spect = int_spect(:,indexL:indexR);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Index exceeds matrix dimensions.
    mean_int = mean(temp_spect);
    threshold = mean_int + 0.1;
    %disp(n);disp(threshold)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if int1 < threshold
        eliminate_pks = [eliminate_pks,n];
    end
end
locs(eliminate_pks) = [];
number_of_peaks = length(locs);

%look at the binary image called bounds at the indices where the intensity
%peaks are to create new spectrums
pk_spectrum = cell(1,length(locs));
for n = 1:length(locs)
    pk_spectrum{n} = Ibounds(:,locs(n));
end

%find the largest sets of straight ones in the peak spectrums and declare
%those indices as the middle region, and then kill off any fringes that
%don't have chosen middle sections
binary_idx = cell(1,length(locs));
[a,b] = size(I);
siz = [a b];
for n = 1:length(locs)
    ones_idx = find(pk_spectrum{n} == 1);
    position_idx = ones(length(ones_idx),1)*locs(n);
    linearInd = sub2ind(siz, ones_idx, position_idx);
    binary_idx{n} = [position_idx, ones_idx, linearInd];
end
fringes_at_locs = zeros(1,length(locs));
for m = 1:length(locs)
    lngth = zeros(1,length(fringe_container));
    for n = 1:length(fringe_container)
        C = intersect(binary_idx{m}(:,3),fringe_container{n});
        temp = length(C);
        lngth(n) = temp;
    end
    [~,max_idx] = max(lngth);
    fringes_at_locs(m) = max_idx;
end
important_fringes_idx = unique(fringes_at_locs);
important_fringes = cell(1,length(important_fringes_idx));
fringe_containers = cell(1,length(important_fringes_idx));
important_fringes_containers = cell(1,length(important_fringes_idx));
for n = 1:length(important_fringes_idx)
    important_fringes_containers{n} = fringe_container{important_fringes_idx(n)};
    important_fringes{n} = fringes_sorted{important_fringes_idx(n)};
    fringe_containers{n} = fringe_container{important_fringes_idx(n)};
end
%figure;plot(important_fringes{1}(:,6))

%visualization
%figure;title(Ipack.fileName);subplot(2,1,1);imshow(I);subplot(2,1,2);imshow(Ibounds2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[x,y]=size(Ibounds);
%X=1:x;
%Y=1:y;
%[xx,yy]=meshgrid(Y,X);
%figure;surf(xx,yy,Ibounds);
%figure;surf(xx,yy,Ibounds2);

%create new cell Ibounds new to store all the pixel points within the area
%of the fringes (from max curve edge to max curve edge)
%Ibounds_new = cell(1,length(important_fringes));
Ifit_data = cell(1,length(important_fringes));
%Ifit_data = [row_index, middle_column_index, average_column_index]
for n = 1:length(important_fringes)
    [c,~] = size(important_fringes{n});
    for m = 1:c
        span = important_fringes{n}(m,4):important_fringes{n}(m,5);
        span = span';
        %temp1 = ones(length(span),1)*important_fringes{n}(m,17);
        %Ibounds_new{n} = [Ibounds_new{n}; temp1, span];
        span_avg = mean(span);
        Ifit_data{n} = [Ifit_data{n};span_avg];
    end
    Ifit_data{n} = [important_fringes{n}(:,17),important_fringes{n}(:,6),Ifit_data{n}];
end

%make two smoothed fits: 1 using the fringe middles and one using the
%fringe average value
Ifit = cell(1,length(important_fringes));
fringe_center = [];
%Ifit = [row_index, middle_column_index, average_column_index]
for n = 1:length(important_fringes)
    temp1 = min(Ifit_data{n}(:,1)):max(Ifit_data{n}(:,1));
    x1 = temp1';
    Ifit{n} = x1;
    x = Ifit_data{n}(:,1);
    y = Ifit_data{n}(:,2);
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 7.2109023599527e-06;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    y1 = fitresult(x1);
    y1 = sgolayfilt(y1,1,31);
    Ifit{n} = [Ifit{n}, y1];
    y = Ifit_data{n}(:,3);
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 7.2109023599527e-06;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    y1 = fitresult(x1);
    y1 = sgolayfilt(y1,1,31);
    Ifit{n} = [Ifit{n}, y1];
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'poly2' );
    [fitresult, gof] = fit( xData, yData, ft );
    y1 = fitresult(x1);
    [~,idx] = min(y1);
    fringe_center = [fringe_center,idx];
    fringe_center = mean(fringe_center);
end

%find the largest continuous section of low slope for each fringe with both
%middle and average fits and then take whichever one is bigger and call
%that the middle, then find the endpoints of the middle spans of all
%relevant fringes and average them to find the indices of the middle
slope_filter = 0.25;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%see if I can change this to compensate for non-flat fringes
span = cell(1,2);
middle = cell(1,length(important_fringes));
mid_bounds = zeros(length(important_fringes),2);
for n = 1:length(important_fringes)
    x = Ifit{n}(:,1);
    y = Ifit{n}(:,2);
    y_slope = gradient(y);
    y_slope = sgolayfilt(y_slope,1,21);
    Y_slope_binary = y_slope;
    Y_slope_binary(abs(y_slope) > slope_filter) = 0;
    Y_slope_binary(Y_slope_binary ~= 0) = 1;
    Y_slope_binary = Y_slope_binary';
    measurements = regionprops(logical(Y_slope_binary), 'Area','PixelIdxList');
    [~,idx] = max([measurements.Area]);
    temp1 = measurements(idx).PixelIdxList;
    span{1} = x(temp1);
    y = Ifit{n}(:,3);
    y_slope = gradient(y);
    y_slope = sgolayfilt(y_slope,1,21);
    Y_slope_binary = y_slope;
    Y_slope_binary(abs(y_slope) > slope_filter) = 0;
    Y_slope_binary(Y_slope_binary ~= 0) = 1;
    Y_slope_binary = Y_slope_binary';
    measurements = regionprops(logical(Y_slope_binary), 'Area','PixelIdxList');
    [~,idx] = max([measurements.Area]);
    temp1 = measurements(idx).PixelIdxList;
    span{2} = x(temp1);
    [~,idx] = max([length(span{1}),length(span{2})]);
    middle{n} = span{idx};
    mid_bounds(n,1) = span{idx}(1);
    mid_bounds(n,2) = span{idx}(length(span{idx}));
end
middleLower = round(max(mid_bounds(:,1)));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed from mean to max and min
middleUpper = round(min(mid_bounds(:,2)));
center = round(mean([middleLower middleUpper]));
%assignin('base','mid_bounds',mid_bounds)

%find average center position
approx_centers = zeros(1,length(important_fringes));
for n = 1:length(important_fringes)
    fa = important_fringes{n};
    fa17 = fa(:,17);
    [temp1,~] = find(fa17 < middleLower);
    [temp2,~] = find(fa17 > middleUpper);
    temp = [temp1;temp2];
    fa(temp,:) = [];
    avg_cent = mean(fa(:,6));
    approx_centers(n) = round(avg_cent);
end

%choose a side length for fringes    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed side length   !!!!!!!!!!!!!
side_length = 75;

%Isolate the down side
slope_filter = -0.35;
%span = cell(1,2);
down = cell(1,length(important_fringes));
down_bounds = zeros(length(important_fringes),2);
%down_bounds = [upper-most row, lower-most row], nx2 array, 1 row per fringe
slope = zeros(length(important_fringes),2);
for n = 1:length(important_fringes)
    x = Ifit{n}(:,1);
    y = Ifit{n}(:,2);
    y_slope = gradient(y);
    y_slope = sgolayfilt(y_slope,1,21);
    Y_slope_binary = y_slope;
    Y_slope_binary(y_slope > slope_filter) = 0;
    Y_slope_binary(Y_slope_binary ~= 0) = 1;
    Y_slope_binary = Y_slope_binary';
    measurements = regionprops(logical(Y_slope_binary), 'Area','PixelIdxList');
    [~,idx] = max([measurements.Area]);
    temp1 = measurements(idx).PixelIdxList;
    down{n} = x(temp1);
    down_bounds(n,2) = down{n}(length(down{n}));
    if (down{n}(length(down{n})) - side_length) > down{n}(1)
        down_bounds(n,1) = down{n}(length(down{n})) - side_length + 1;
        temp1 = temp1(end-(side_length-1):end);
    else
        down_bounds(n,1) = down{n}(1);
    end
    x = x(temp1);
    down{n} = x;
    y = y(temp1);
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    y1 = fitresult(x);
    slope(n,1) = ((y1(length(y1))-y1(1))/(x(length(x))-x(1)))*-1;
end

%Isolate the up side
slope_filter = 0.35;
%span = cell(1,2);
up = cell(1,length(important_fringes));
up_bounds = zeros(length(important_fringes),2);
%up_bounds = [lower-most row, upper-most row], nx2 array, 1 row per fringe
for n = 1:length(important_fringes)
    x = Ifit{n}(:,1);
    y = Ifit{n}(:,2);
    y_slope = gradient(y);
    y_slope = sgolayfilt(y_slope,1,21);
    Y_slope_binary = y_slope;
    Y_slope_binary(y_slope < slope_filter) = 0;
    Y_slope_binary(Y_slope_binary ~= 0) = 1;
    Y_slope_binary = Y_slope_binary';
    measurements = regionprops(logical(Y_slope_binary), 'Area','PixelIdxList');
    [~,idx] = max([measurements.Area]);
    temp1 = measurements(idx).PixelIdxList;
    up{n} = x(temp1);
    up_bounds(n,1) = up{n}(1);
    if (up{n}(1) + side_length) < up{n}(length(up{n}))
        up_bounds(n,2) = up{n}(1) + side_length - 1;
        temp1 = temp1(1:side_length);
    else
        up_bounds(n,2) = up{n}(length(up{n}));
    end
    x = x(temp1);
    up{n} = x;
    y = y(temp1);
    [xData, yData] = prepareCurveData( x, y );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    y1 = fitresult(x);
    slope(n,2) = ((y1(length(y1))-y1(1))/(x(length(x))-x(1)));
end

%determine which fringes are odd/even, store in vector odd or even
[a,~] = size(slope);
slope1 = zeros(a,1);
for n = 1:a
    slope1(n) = mean([slope(n,1), slope(n,2)]);
end
odd_or_even = zeros(a,1);   %odds are 1, evens are 2
if mod(a,2) == 0
    k = a/2;
    [~,idx] = maxk(slope1,k);
    odd_or_even(idx) = 1;
    odd_or_even(odd_or_even ~= 1) = 2;
else
    k = (a/2)-0.5;
    [~,idx1] = maxk(slope1,k);
    [~,idx2] = mink(slope1,k);
    odd_or_even(idx1) = 1;
    odd_or_even(idx2) = 2;
    idx = find(odd_or_even == 0);
    temp = slope1(idx);
    temp1 = slope1(idx1);
    temp1 = min(temp1);
    temp2 = slope1(idx2);
    temp2 = min(temp2);
    b = abs(temp-temp1);
    c = abs(temp-temp2);
    bc = [b c];
    [~,bc] = min(bc);
    if bc == 1
        odd_or_even(idx) = 1; %odd
    else
        odd_or_even(idx) = 2; %even
    end
end
slope = [slope1, odd_or_even];
%mod(N,2)==0 % iseven
%mod(N,2)==1 % isodd
%fix(N)==N % isint

%find any gamma fringes or shadows, store positions in double fringes
double_fringes = [];
[a,~] = size(slope);
for n=1:(a-1)
    if slope(n,2) == slope(n+1,2)
        double_fringes = [double_fringes, (n+1)];
    end
end
%info for Ipack: important fringes, approximate centers, middle points,
%Ifit, sides,odd/even labels
%
%inform user if any shadow fringes are found
if length(double_fringes) > 0
    disp('making an odd - even deletion!!! fringes numbers to be eliminated:')
    disp(double_fringes)
end
double_fringes = [];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%eliminate all the fringes with indices of the double_fringes array
slope(double_fringes,:) = [];   %[average edge slope, odd/even]
up_bounds(double_fringes,:) = [];   %[first/last x indices of upward edges]
down_bounds(double_fringes,:) = []; %[first/last x indices of downward edges]
approx_centers = approx_centers';   %center of fringe in terms of the wavelength axis
approx_centers(double_fringes,:) = [];
up(double_fringes) = [];
down(double_fringes) = [];
middle(double_fringes) = [];
Ifit(double_fringes) = [];
important_fringes(double_fringes) = [];
fringe_containers(double_fringes) = [];
%for n = 1:length(double_fringes)
%    up(double_fringes(n)) = [];
%    down(double_fringes(n)) = [];
%    middle(double_fringes(n)) = [];
%    Ifit(double_fringes(n)) = [];
%    important_fringes(double_fringes(n)) = [];
%    fringe_containers(double_fringes(n)) = [];
%end
%other important info: middleLower, middleUpper, center, lower_bound,
%upper_bound

%figure out where to crop the image for aligning the middles and where to
%crop the image for the error between odd/even meas
%
%use middleLower:middleUpper for cropping middle and the minimum of 2nd
%column of up bounds:max of down bounds first column for cropping with
%edges
lower_bound = max(down_bounds(:,1));
upper_bound = min(up_bounds(:,2));
%assignin('base','up_bounds',up_bounds)
%assignin('base','down_bounds',down_bounds)
%
%shorten fringe array to align with crops and change upper/lower bound if
%need be
fringe_array = important_fringes;
fringe_array_sides = cell(1,length(fringe_array));
fringe_array_middle = cell(1,length(fringe_array));
%disp(upper_bound)
for n = 1:length(fringe_array)
    %for fringe array middle
    [~,idx,~] = intersect(fringe_array{n}(:,17),[middleLower, middleUpper]);
    fringe_array_middle{n} = fringe_array{n}(idx(1):idx(2),:);
    %
    %for fringe_array sides
    [~,idx,~] = intersect(fringe_array{n}(:,17),lower_bound);
    if length(idx) == 1
        %idx1 = fringe_array{n}(idx,17);
        %do nothing
    elseif length(idx) == 0
        if min(fringe_array{n}(:,17)) > lower_bound
            idx1 = min(fringe_array{n}(:,17));
            lower_bound = idx1;
        elseif min(fringe_array{n}(:,17)) < lower_bound
            temp = fringe_array{n}(:,17);
            temp(find(temp < lower_bound)) = [];
            idx1 = min(temp);
            lower_bound = idx1;
        end
    end
    [~,idx,~] = intersect(fringe_array{n}(:,17),upper_bound);
    if length(idx) == 1
        %idx2 = fringe_array{n}(idx,17);
        %do nothing
    elseif length(idx) == 0
        if max(fringe_array{n}(:,17)) < upper_bound
            idx2 = max(fringe_array{n}(:,17));
            upper_bound = idx2;
        elseif max(fringe_array{n}(:,17)) > upper_bound
            temp = fringe_array{n}(:,17);
            temp(find(temp > upper_bound)) = [];
            idx2 = max(temp);
            upper_bound = idx2;
        end
    end
end
%disp([lower_bound 1])
%disp([upper_bound 2])
%assignin('base','fringe_array',fringe_array)
for n = 1:length(fringe_array)
    [~,idx,~] = intersect(fringe_array{n}(:,17),[upper_bound, lower_bound]);
    %disp(n)
    %disp(idx)
    fringe_array_sides{n} = fringe_array{n}(idx(1):idx(2),:);
end
%
%crop the images using previously defined maximum and minimum bounds
Icrop_orig = Ipack.crop_smoothed;
Icrop_middle = Icrop_orig(middleLower:middleUpper,:);
Icrop_sides = Icrop_orig(lower_bound:upper_bound,:);
%
%record crop indices for middle and sides crops. These will tell you what
%row in the original image the crop starts at so that to translate the
%index from the croppped to the uncropped image you just need to add the
%crop index to the cropped image index.
middle_crop_index = middleLower;
sides_crop_index = lower_bound;
%
%find the edges of the fringe containers
[a,b] = size(Ibounds2);
siz = [a b];
fringe_container2 = cell(1,length(fringe_containers));
%fringe_container2: [row, lower column, upper column]
for n = 1:length(fringe_containers)
    [I1,J1] = ind2sub(siz,fringe_containers{n});
    idx = [I1, J1];
    idx = sortrows(idx);
    idx(idx(:,1) > upper_bound,:) = [];
    idx(idx(:,1) < lower_bound,:) = [];
    min_idx = min(idx(:,1));
    max_idx = max(idx(:,1));
    span = min_idx:max_idx;
    container_edges = zeros(length(span),3);
    for l = 1:length(span)
        m = span(l);
        idx2 = find(idx(:,1) == m);
        temp1 = min(idx(idx2,2));
        temp2 = max(idx(idx2,2));
        container_edges(l,:) = [m temp1 temp2];
    end
    fringe_container2{n} = container_edges;
end
%
%recalculate indices to work with enhanced, cropped image
%
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
%k(14) if more than 2 ridges, the indices of the 1st of the two highest ridges - i.e. refers to fringe array not to image
%l(15) if more than 2 ridges, the indices of the 2nd of the two highest ridges
%m(16) number of ridges within the fringe
%n(17) row fringe is situated in
%o(18) linear index of fringe middle
%
%change all old indices to enhanced 1/10th indices
%where f(x) is new enhanced and cropped indices
%for horizontal indices: f(x) = (x*10)-9
%for vertical indices: f(x) = [x-(cropping_index-1)]*10 - 9
fa_horiz_idx = [1 2 4 5 6 7 8 9 10 11 12 13];
for n = 1:length(fringe_array)
    fringe_array_sides{n}(:,fa_horiz_idx) = (fringe_array_sides{n}(:,fa_horiz_idx).*10)-9;
    [a,b] = find(fringe_array_sides{n} == -9);
    fringe_array_sides{n}(a,b) = 0;
    fringe_array_sides{n}(:,17) = (fringe_array_sides{n}(:,17) - (sides_crop_index - 1)).*10 - 9;
    fringe_array_middle{n}(:,fa_horiz_idx) = (fringe_array_middle{n}(:,fa_horiz_idx).*10)-9;
    [a,b] = find(fringe_array_middle{n} == -9);
    fringe_array_middle{n}(a,b) = 0;
    fringe_array_middle{n}(:,17) = (fringe_array_middle{n}(:,17) - (middle_crop_index - 1)).*10 - 9;
    approx_centers_enhanced = (approx_centers.*10)-9;
    fringe_container2{n}(:,2:3) = (fringe_container2{n}(:,2:3) .* 10) - 9;
    fringe_container2{n}(:,1) = (fringe_container2{n}(:,1) - (sides_crop_index - 1)).*10 - 9;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%assignin('base','approx_centers',approx_centers)
%
%enhance the cropped images to 0.1 the resolution
Icrop_middle = enhanceImage(Icrop_middle);
Icrop_sides = enhanceImage(Icrop_sides);
%
%make an enhanced fit of the fringe container edges to match the enhanced
%resolution images
fringe_containers = cell(1,length(fringe_container2));
for n = 1:length(fringe_container2)
    x = fringe_container2{n}(:,1);
    y = fringe_container2{n}(:,2);
    z = fringe_container2{n}(:,3);
    x1 = min(x):max(x);
    x1 = x1';
    [xData, yData] = prepareCurveData( x, y );
    ft = 'linearinterp';
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    y1 = round(fitresult(x1));
    [xData, zData] = prepareCurveData( x, z );
    ft = 'linearinterp';
    [fitresult, gof] = fit( xData, zData, ft, 'Normalize', 'on' );
    z1 = round(fitresult(x1));
    fringe_containers{n} = [x1 y1 z1];
end

%calculate highest and lowest points for each fringe container
%fringe_containers: [row, lower column, upper column]
fringe_container_bounds = zeros(length(fringe_containers),2);
%fringe_container_bounds: [lowest point of fringe container n, highest point of fringe continaer n]
for n = 1:length(fringe_containers)
    fc = fringe_containers{n};
    fringe_container_bounds(n,1) = min(fc(:,2));
    fringe_container_bounds(n,2) = max(fc(:,3));
end

%save information in Ipack
%Ipack.fringe_fit = Ifit;
%Ipack.up = up;
%Ipack.down = down;
%Ipack.middle = middle;
Ipack.fringe_array = important_fringes;
Ipack.middle_crop_index = middle_crop_index;
Ipack.sides_crop_index = sides_crop_index;
Ipack.crop_middle = Icrop_middle;
Ipack.crop_sides = Icrop_sides;
Ipack.approx_centers = approx_centers_enhanced;
Ipack.fringe_array_sides = fringe_array_sides;
Ipack.fringe_array_middle = fringe_array_middle;
Ipack.fringe_containers = fringe_containers;
Ipack.fringe_container_bounds = fringe_container_bounds;


%show pictures of how successful everything was
%visualization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%n = 2;
%x = Ifit{n}(:,1);
%y = Ifit{n}(:,2);
%y_slope = gradient(y);
%y_slope = sgolayfilt(y_slope,1,21);
%y_curv = gradient(y_slope);
%figure;plot(x,y)
%figure;plot(x,y_slope);
%figure;plot(x,y_curv);
%
%figure;imshow(Icrop_middle)
%figure;imshow(Icrop_sides)
%
figure;subplot(2,1,1);imshow(I);title(Ipack.fileName);
hold on
for n = 1:length(Ifit)
    fit1 = Ifit{n};
    down1 = down{n};
    middle1 = middle{n};
    up1 = up{n};
    [~,idx_down,~] = intersect(fit1(:,1),down1);
    [~,idx_middle,~] = intersect(fit1(:,1),middle1);
    [~,idx_up,~] = intersect(fit1(:,1),up1);
    down2 = fit1(idx_down,1:2);
    middle2 = fit1(idx_middle,1:2);
    up2 = fit1(idx_up,1:2);
    line(down2(:,2),down2(:,1),'color','r');
    line(middle2(:,2),middle2(:,1),'color','b');
    line(up2(:,2),up2(:,1),'color','r');
end
hold off
num_pks = num2str(number_of_peaks);
plot_tit = strcat('#peaks:',num_pks);
subplot(2,1,2);imshow(Ibounds2);title(plot_tit);

%assignin('base','lower_bound',lower_bound);
%assignin('base','upper_bound',upper_bound);
%assignin('base','fringe_container',fringe_container);
%assignin('base','Ibounds2',Ibounds2);
%Ipack.lb = lower_bound;
%Ipack.ub = upper_bound;
%Ipack.fringe_container = fringe_containers;
%Ipack.Ibinary = Ibounds2;

end