%compare fringe areas of same shape as 1st fringe and analyze the height
%difference error within them

%zero_curv_bounds_smoothed, 1xn cell, [row, low column, high column]
%
%comparison_array, (n-1)x6 array
%comparison array: [(1)fringe n idx,(2)fringe n+1 idx,(3)fringe n center,(4)fringe n+1 new center,(5)alignment error,(6)fringe n+1 offset from original center]
%2nd fringes original center = (new center) - (offset)
%
%I = Ipack.crop_sides
%
%height_adjust is height adjusted during fringe alignment, nx1 vector for n
%fringes

%expand zero_curv_bounds_smoothed to a matrix that includes bounds for both
%fringes n and n+1, called bounds_comparison
[a,~] = size(comparison_array);
bounds_comparison = cell(1,a);
for n = 1:a
    curv0 = zero_curv_bounds_smoothed{n};
    curv0(:,1) = [];
    cent_1 = comparison_array(n,3);
    cent_2 = comparison_array(n,4);
    curv2 = zeros(length(curv0),2);
    for m = 1:length(curv0)
        bound1 = curv0(m,1);
        bound2 = curv0(m,2);
        temp1 = bound1 - cent_1;
        temp2 = bound2 - cent_1;
        bound1_new = cent_2 + temp1;
        bound2_new = cent_2 + temp2;
        curv2(m,:) = [bound1_new bound2_new];
    end
    bounds_comparison{n} = [curv0 curv2];
end
%n = 1;
%x = 1:length(bounds_comparison{n});
%x = x';
%figure;plot(x,bounds_comparison{n}(:,1),x,bounds_comparison{n}(:,2),x,bounds_comparison{n}(:,3),x,bounds_comparison{n}(:,4))


error_plots = cell(1,length(bounds_comparison));
for n = 1:length(bounds_comparison)
    bounds = bounds_comparison{n};
    err_dat = zeros(length(bounds),1);
    for m = 1:length(bounds)
        span1 = I(m,bounds(m,1):bounds(m,2));
        span2 = I(m,bounds(m,3):bounds(m,4));
        span1 = span1 + height_adjust(n);
        span2 = span2 + height_adjust(n+1);
        span_error = span1 - span2;
        span_error = span_error.^2;
        num_elements = length(span_error);
        span_error = sum(span_error);
        span_error = span_error / num_elements;
        span_error = span_error^0.5;
        err_dat(m) = span_error;
    end
    error_plots{n} = err_dat;
    figure;plot(err_dat)
end
