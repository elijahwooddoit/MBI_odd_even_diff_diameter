%move the fringes within the zero curv bounds to the blank image shifted
%up corresponding to the fringe alignment

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

%fill in two binary images with shifted image within bounds, shifted to the
%correct height
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
        span = span + height_adjust(n);
        span_shift = bound_right - bound_left;
        temp = cent1 - bound_left;
        new_start = cent2 - temp;
        IB1(m,new_start:(new_start+span_shift)) = span;
        %for IB2
        bound_left = bounds_comparison{n}(m,3);
        bound_right = bounds_comparison{n}(m,4);
        span = I(m,bound_left:bound_right);
        span = span + height_adjust(n+1);
        span_shift = bound_right - bound_left;
        temp = mid1 - bound_left;
        new_start = mid2 - temp;
        IB2(m,new_start:(new_start+span_shift)) = span;
    end
    comparison_images{n}{1} = IB1;
    comparison_images{n}{2} = IB2;
    %[x,y]=size(Ibinary_blank);
    %X=1:x;
    %Y=1:y;
    %[xx,yy]=meshgrid(Y,X);
    %figure;surf(xx,yy,IB1,'linestyle','none');
    %figure;surf(xx,yy,IB2,'linestyle','none');
end

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
    figure;plot(err_dat)
end