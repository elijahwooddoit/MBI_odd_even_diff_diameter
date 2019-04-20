function Ipack = ridge_find(Ipack)

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

%find peaks in the intensity signal and make sure they're tall enough
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
important_fringes_containers = cell(1,length(important_fringes_idx));
for n = 1:length(important_fringes_idx)
    important_fringes_containers{n} = fringe_container{important_fringes_idx(n)};
    important_fringes{n} = fringes_sorted{important_fringes_idx(n)};
end
%figure;plot(important_fringes{1}(:,6))

%visualization
figure;title(Ipack.fileName);subplot(2,1,1);imshow(I);subplot(2,1,2);imshow(Ibounds2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
middleL = round(mean(mid_bounds(:,1)));
middleR = round(mean(mid_bounds(:,2)));
center = round(mean([middleL middleR]));

%visualization
%n = 2;
%x = Ifit{n}(:,1);
%y = Ifit{n}(:,2);
%y_slope = gradient(y);
%y_slope = sgolayfilt(y_slope,1,21);
%y_curv = gradient(y_slope);
%figure;plot(x,y)
%figure;plot(x,y_slope);
%figure;plot(x,y_curv);