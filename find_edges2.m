function Ipack = find_edges2(Ipack)

%make rectangular images of fringe middles to find center, and then
%immediately find the center of each fringe
approx_centers = Ipack.approx_centers;
fringe_array_middle = Ipack.fringe_array_middle;
I = Ipack.crop_middle;
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

%determine the maxiumum distances from the center to each edge of the
%fringe, then find the longest distances among all the fringes and keep
%those
fringe_width = zeros(length(fringe_array_middle),2);
for n = 1:length(fringe_array_middle)
    min_bound = min(fringe_array_middle{n}(:,7));%%%%%%%%% changed fringe_array indices from 4/5 to 7/8 
    max_bound = max(fringe_array_middle{n}(:,8));
    cent = approx_centers(n);
    fringe_width(n,1) = cent - min_bound;
    fringe_width(n,2) = max_bound - cent;
end
final_width = [max(fringe_width(:,1)), max(fringe_width(:,2))];

%using the width from above, determine what indices to crop each fringe at
%and where their centers will be in the new crops
crop_indices = zeros(length(fringe_array_middle),2);
crop_centers = zeros(1,length(fringe_array_middle));
[~,b] = size(I);
for n = 1:length(fringe_array_middle)
    if approx_centers(n)-final_width(1) > 0
        crop_indices(n,1) = approx_centers(n)-final_width(1);
        crop_centers(n) = final_width(1) + 1;
    else
        crop_indices(n,1) = 1;
        crop_centers(n) = final_width(1) + (approx_centers(n)-final_width(1));
    end
    if approx_centers(n)+final_width(2) < b
        crop_indices(n,2) = approx_centers(n)+final_width(2);
    else
        crop_indices(n,2) = b;
    end
end

%shift all fringe maximums up to 1 an save the shift amount
cropped_middles = cell(1,length(fringe_array_middle));
height_adjust = zeros(length(fringe_array_middle),1);
for n = 1:length(fringe_array_middle)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changing this from 0-1 to shift up to 1
    cropped_middles{n} = I(:,crop_indices(n,1):crop_indices(n,2));
    %temp1 = max(max(cropped_middles{n}));
    temp1 = mean(mean(cropped_middles{n}));
    temp2 = 1 - temp1;
    cropped_middles{n} = cropped_middles{n} + temp2;
    height_adjust(n) = temp2;
    %[x,y]=size(cropped_middles{n});
    %X=1:x;
    %Y=1:y;
    %[xx,yy]=meshgrid(Y,X);
    %figure;surf(xx,yy,cropped_middles{n},'linestyle','none');title(num2str(n));
end
%
%deleted version:
%normalize all the fringes to an intensity range of 0-1 and save the
%normalization info for each fringe for use later on
%normalization = zeros(length(fringe_array_middle),2);
%normalization: [amount to subtract, amount to multiply by]
%cropped_middles = cell(1,length(fringe_array_middle));
%for n = 1:length(fringe_array_middle)
%    cropped_middles{n} = I(:,crop_indices(n,1):crop_indices(n,2));
%    temp1 = min(min(cropped_middles{n}));
%    cropped_middles{n} = cropped_middles{n} - temp1;
%    temp2 = 1/max(max(cropped_middles{n}));
%    cropped_middles{n} = cropped_middles{n} * temp2;
%    normalization(n,1) = temp1;
%    normalization(n,2) = temp2;
%end

%analyze error between all sequential fringes to determine which sequential
%set has the least amount of error, and which fringes to use as base,
%base+1
comparison_array = zeros(length(fringe_array_middle) - 1,6);
%comparison array: [(1)fringe n idx,(2)fringe n+1 idx,(3)fringe n center,(4)fringe n+1 new center,(5)alignment error,(6)fringe n+1 offset from original center]
%2nd fringes original center = (new center) - (offset)
for n = 1:(length(fringe_array_middle) - 1)
    Im1 = cropped_middles{n};
    Im2 = cropped_middles{n+1};
    center1 = crop_centers(n);
    center2 = crop_centers(n+1);
    cent_err_offset = fringe_aligner(Im1,center1,Im2,center2);
    comparison_array(n,1) = n;
    comparison_array(n,2) = n+1;
    comparison_array(n,3) = approx_centers(n);
    comparison_array(n,4) = approx_centers(n+1) + cent_err_offset(3);
    comparison_array(n,5:6) = cent_err_offset(2:3);
end

%mod(N,2)==0 % iseven
%mod(N,2)==1 % isodd
%fix(N)==N % isint
clear Im1 Im2 I fringe_array_middle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%isolate each fringe with the sides included and build error plots of each
%consecutive pair of fringes

I = Ipack.crop_sides;
fringe_containers = Ipack.fringe_containers;
fringe_array = Ipack.fringe_array_sides;
fringe_container_bounds = Ipack.fringe_container_bounds;

%fringe_container_bounds: [(horizontally speaking) lowest point of fringe container n, highest point of fringe continaer n]
%
%fringe_containers: [row, lower column, upper column]
%
%%comparison array: [(1)fringe n idx,(2)fringe n+1 idx,(3)fringe n center,(4)fringe n+1 new center,(5)alignment error,(6)fringe n+1 offset from original center]
%2nd fringes original center = (new center) - (offset)
%
%approx_centers: vertical vector of original center approximations
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

%create a blank binary image that is big enough to fit any of the fringes
%within it, with its center as 'center_position'
dist_from_centers = zeros(length(fringe_array),2);
shift_left2 = min(comparison_array(find(comparison_array(:,6) < 0),6));
if isempty(shift_left2) == 1
    shift_left2 = 0;
end
shift_right2 = max(comparison_array(find(comparison_array(:,6) > 0),6));
if isempty(shift_right2) == 1
    shift_right2 = 0;
end
for n = 1:length(fringe_array)
    dist_from_centers(n,1) = approx_centers(n) - fringe_container_bounds(n,1);
    dist_from_centers(n,2) = fringe_container_bounds(n,2) - approx_centers(n);
end
shift_left1 = max(dist_from_centers(:,1));
shift_right1 = max(dist_from_centers(:,2));
shift_left = shift_left1 + shift_left2 + 50;
shift_right = shift_right1 + shift_right2 + 50;
[a,~] = size(I);
Ibinary_blank = zeros(a,shift_left + shift_right + 1);
center_position = shift_left + 1;

%find the zero curvature edges for the new image
zero_curv_bounds = cell(1,length(fringe_array));
zero_curv_bounds_fit = cell(1,length(fringe_array));
zero_curv_bounds_smoothed = cell(1,length(fringe_array));
%zero_curv_bounds_smoothed, 1xn cell, [row, low column, high column]
for n = 1:length(fringe_array)
    x = fringe_array{n}(:,17);
    y = fringe_array{n}(:,7);
    [xData, yData] = prepareCurveData( x, y );
    ft = 'linearinterp';
    [fitresult, ~] = fit( xData, yData, ft, 'Normalize', 'on' );
    x1 = min(x):max(x);
    x1= x1';
    y1 = fitresult(x1);
    %figure;plot(x1,y1)
    z = fringe_array{n}(:,8);
    [xData, zData] = prepareCurveData( x, z );
    ft = 'linearinterp';
    [fitresult, gof] = fit( xData, zData, ft, 'Normalize', 'on' );
    x1 = min(x):max(x);
    x1= x1';
    z1 = fitresult(x1);
    %figure;plot(x1,y1)
    zero_curv_bounds_fit{n} = [x1 y1 z1];
    [a,b] = size(zero_curv_bounds_fit{n});
    zero_curv_bounds{n} = zeros(a,b);
    for m = 1:length(x1)
        zero_curv_bounds{n}(m,1) = m;
        y = I(fringe_containers{n}(m,1),fringe_containers{n}(m,2):fringe_containers{n}(m,3));
        y = sgolayfilt(y,2,21);
        y_deriv1 = gradient(y);
        y_deriv1 = sgolayfilt(y_deriv1,2,21);
        y_deriv2 = gradient(y_deriv1);
        y_deriv2 = sgolayfilt(y_deriv2,2,21);
        x = 1:length(y);
        x = x';
        z = zeros(length(y),1);
        [~,~,iout,~] = intersections(x,y_deriv2,x,z);%x,y_deriv2,x,z
        [~,idx] = min(abs(iout-(zero_curv_bounds_fit{n}(m,2)-fringe_containers{n}(m,2)+1)));
        minVal1 = round(iout(idx));
        [~,idx] = min(abs(iout-(zero_curv_bounds_fit{n}(m,3)-fringe_containers{n}(m,2)+1)));
        minVal2 = round(iout(idx));
        val1 = minVal1 + fringe_containers{n}(m,2) - 1;
        val2 = minVal2 + fringe_containers{n}(m,2) - 1;
        zero_curv_bounds{n}(m,2:3) = [val1 val2];
    end
    %original commands:
    %temp1 = round(sgolayfilt(zero_curv_bounds{n}(:,2),3,151));
    %temp2 = round(sgolayfilt(zero_curv_bounds{n}(:,3),3,151));
    %zero_curv_bounds_smoothed{n} = [zero_curv_bounds{n}(:,1) temp1 temp2];
    
    %new commands playing around:
    poly_deg = 1;
    filt_window = 301;
    tt1 = num2str(poly_deg);
    tt2 = num2str(filt_window);
    plot_tit = strcat(tt1, '%', tt2);
    temp1 = round(sgolayfilt(zero_curv_bounds{n}(:,2),poly_deg,filt_window));
    temp2 = round(sgolayfilt(zero_curv_bounds{n}(:,3),poly_deg,filt_window));
    zero_curv_bounds_smoothed{n} = [zero_curv_bounds{n}(:,1) temp1 temp2];
    %figure;plot(zero_curv_bounds_smoothed{n}(:,1),zero_curv_bounds_smoothed{n}(:,2),zero_curv_bounds_smoothed{n}(:,1),zero_curv_bounds_smoothed{n}(:,3),...
    %    zero_curv_bounds{n}(:,1),zero_curv_bounds{n}(:,2),zero_curv_bounds{n}(:,1),zero_curv_bounds{n}(:,3));title(plot_tit)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the fringe ridges within the zero_curv_bounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make smoothed binary image of the fringe within zero curv bounds for
%calculating the error
%
%zero_curv_bounds_smoothed, 1xn cell, [row, low column, high column]
%
%comparison_array, (n-1)x6 array
%comparison array: [(1)fringe n idx,(2)fringe n+1 idx,(3)fringe n center,(4)fringe n+1 new center,(5)alignment error,(6)fringe n+1 offset from original center]
%2nd fringes original center = (new center) - (offset)
%
%I = Ipack.crop_sides
%
%Ibinary_blank = zeros(a,shift_left + shift_right + 1);
%center_position = shift_left + 1;
%
%height_adjust is height adjusted during fringe alignment, nx1 vector for n
%fringes

%expand zero_curv_bounds_smoothed to a matrix that includes bounds for both
%fringes n and n+1, called bounds_comparison
[a,~] = size(comparison_array);
bounds_comparison = cell(1,a);
%bounds_comparison: [n fringe left bound, n fringe right bound, n+1 fringe left bound, n+1 fringe right bound]
for n = 1:a
    curv0 = zero_curv_bounds_smoothed{n};
    curv0(:,1) = [];
    curv2 = zero_curv_bounds_smoothed{n+1};
    curv2(:,1) = [];
    bounds_comparison{n} = [curv0 curv2];
end
%n = 1;
%x = 1:length(bounds_comparison{n});
%x = x';
%figure;plot(x,bounds_comparison{n}(:,1),x,bounds_comparison{n}(:,2),x,bounds_comparison{n}(:,3),x,bounds_comparison{n}(:,4))

%fill in two binary images with shifted image within bounds, all points
%within bounds have intensity 1, everything else has intensity 0
[a,~] = size(comparison_array);
comparison_images = cell(1,a);
%comparison_array: cell length of comparison array, where each cell is a
%subdivided into two sub-cells containing images of (n) and (n+1) fringes
for n = 1:a
    comparison_images{n} = cell(1,2);
    IB1 = Ibinary_blank;
    IB2 = Ibinary_blank;
    %for IB1
    cent1 = comparison_array(n,3);
    cent2 = center_position;
    %for IB2
    mid1 = comparison_array(n,4);
    mid2 = center_position;
    for m = 1:length(bounds_comparison{n})
        %for IB1
        bound_left = bounds_comparison{n}(m,1);
        bound_right = bounds_comparison{n}(m,2);
        span = I(m,bound_left:bound_right);
        span(:,:) = 1;
        %span = span + height_adjust(n);
        span_shift = bound_right - bound_left;
        temp = cent1 - bound_left;
        new_start = cent2 - temp;
        IB1(m,new_start:(new_start+span_shift)) = span;
        %for IB2
        bound_left = bounds_comparison{n}(m,3);
        bound_right = bounds_comparison{n}(m,4);
        span = I(m,bound_left:bound_right);
        span(:,:) = 1;
        %span = span + height_adjust(n+1);
        span_shift = bound_right - bound_left;
        temp = mid1 - bound_left;
        new_start = mid2 - temp;
        IB2(m,new_start:(new_start+span_shift)) = span;
    end
    comparison_images{n}{1} = IB1;
    comparison_images{n}{2} = IB2;
    %figure;imshow(IB1,'InitialMagnification','fit')
    %figure;imshow(IB2,'InitialMagnification','fit')
    %[x,y]=size(Ibinary_blank);
    %X=1:x;
    %Y=1:y;
    %[xx,yy]=meshgrid(Y,X);
    %figure;surf(xx,yy,IB1,'linestyle','none');
    %figure;surf(xx,yy,IB2,'linestyle','none');
end

%calculate error by overlaying binary images and summing differences
%between rows
error_plots = cell(1,length(bounds_comparison));
for n = 1:length(comparison_images)
    IB1 = comparison_images{n}{1};
    IB2 = comparison_images{n}{2};
    [b,~] = size(IB1);
    err_dat = zeros(b,1);
    for m = 1:b
        span1 = IB1(m,:);
        span2 = IB2(m,:);
        span_error = span1 - span2;
        span_error = span_error.^2;
        num_elements = length(span_error);
        span_error = sum(span_error);
        span_error = span_error / num_elements;
        span_error = span_error^0.5;
        err_dat(m) = span_error;
    end
    error_plots{n} = err_dat;
    %figure;plot(err_dat)
end

%visualization of error calculations
%
%create a plot of overlayed fringes to get a good visual undertanding of
%how the smoothing is working
[a,~] = size(comparison_array);
shifted_bounds = cell(1,a);
%comparison_array: cell length of comparison array, where each cell is a
%subdivided into two sub-cells containing images of (n) and (n+1) fringes
for n = 1:a
    [b,c] = size(bounds_comparison{n});
    new_bounds = zeros(b,c);
    %for fringe (n)
    cent1 = comparison_array(n,3);
    cent2 = center_position;
    %for fringe (n+1)
    mid1 = comparison_array(n,4);
    mid2 = center_position;
    for m = 1:length(bounds_comparison{n})
        %ffor fringe (n)
        bound_left = bounds_comparison{n}(m,1);
        bound_right = bounds_comparison{n}(m,2);
        span_shift = bound_right - bound_left;
        temp = cent1 - bound_left;
        new_start = cent2 - temp;
        new_end = (new_start+span_shift);
        new_bounds(m,1) = new_start;
        new_bounds(m,2) = new_end;
        %for fringe (n+1)
        bound_left = bounds_comparison{n}(m,3);
        bound_right = bounds_comparison{n}(m,4);
        span_shift = bound_right - bound_left;
        temp = mid1 - bound_left;
        new_start = mid2 - temp;
        new_end = (new_start+span_shift);
        new_bounds(m,3) = new_start;
        new_bounds(m,4) = new_end;
    end
    shifted_bounds{n} = new_bounds;
    %x = 1:b;x = x';
    %figure;plot(x,new_bounds(:,1),'r',x,new_bounds(:,2),'r',x,new_bounds(:,3),'b',x,new_bounds(:,4),'b')
    %figure;plot(error_plots{n})
end

%figure creation and final diameter calculation
[a,~] = size(comparison_array);
name = Ipack.fileName;
for n = 1:a
    [b,~] = size(shifted_bounds{n});
    x = 1:b;x = x';
    temp = num2str(n);
    tit1 = strcat('Outlines; file:',name,'; comp#:',temp);
    %tit2 = strcat('Error; file:',name,'; comp#:',temp);
    error_data = error_plots{n};
    [edges, diameter_odd_even] = errorPlotAnalyzer2(error_data, n, name);
    %figure;plot(x,shifted_bounds{n}(:,1),'r',x,shifted_bounds{n}(:,2),'r',x,shifted_bounds{n}(:,3),'b',x,shifted_bounds{n}(:,4),'b');title(tit1);
    %figure;plot(error_plots{n});title(tit2);
end
    

Ipack.diameter_odd_even = diameter_odd_even;

end


