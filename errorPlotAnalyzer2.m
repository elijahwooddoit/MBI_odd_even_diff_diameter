function [edges, diameter_odd_even] = errorPlotAnalyzer2(error_data, comparison_number, name)

%error_data = error_plots{4};

%normalize error to maximum = 1, minimum = 0, and smooth the error curve,
%and find the gradient
fringe_error = error_data;
x = 1:length(error_data);
y = error_data;
Miny = min(y);
y = y-Miny;
Maxy = max(y);
temp = 1/Maxy;
y = y.*temp;
y_fit = sgolayfilt(y,1,201);
y_deriv = gradient(y_fit);

%find the maximum and minimum slopes of the error curves, which should
%correspond to the edges
[M1,I1] = min(y_deriv);
[M2,I2] = max(y_deriv);

%Find the borders of the down edge
w = 1;
I = I1;
while w == 1    %add in a limiter at I = 1
    M = y_deriv(I);
    if I == 1
        w=0;
    elseif M < -0.001
        I = I+1;
    else
        w = 0;
    end
end
down_max = I;
w = 1;
I = I1;
while w == 1
    M = y_deriv(I);
    if I == 1
        w=0;
    elseif M < -0.001
        I = I-1;
    else
        w = 0;
    end
end
down_min = I;

%find the borders of the up edge
w = 1;
I = I2;
while w == 1
    M = y_deriv(I);
    if M > 0.001
        I = I-1;
    else
        w = 0;
    end
end
up_min = I;
w = 1;
I = I2;
while w == 1
    M = y_deriv(I);%add in a limiter at I equals length of y_deriv
    if I == length(y_deriv)
        w = 0;
    elseif M > 0.001
        I = I+1;
    else
        w = 0;
    end
end
up_max = I;

%define edges
down_errorX = x(down_min:down_max);
down_errorY = y(down_min:down_max);
up_errorX = x(up_min:up_max);
up_errorY = y(up_min:up_max);

%avg height fit of flat part
flat_error = y(down_max:up_min);
mean_error = mean(flat_error);
flat_fit_x = x(down_min:up_max);
flat_fit_y = ones(length(flat_fit_x),1);
flat_fit_y = flat_fit_y*mean_error;

%fit down edge
[xData, yData] = prepareCurveData( down_errorX, down_errorY );
%ft = fittype( 'poly2' );
%[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 6.143039265005923E-7;
[fitresult, gof] = fit( xData, yData, ft, opts );

%find where fit reaches zero
I = down_min;
M = 1;
while M > mean_error
    M = fitresult(x(I));
    I = I+1;
end
down_fit_bot = I;
%make down edge fit line
down_fit_x = x(down_min:down_fit_bot);
down_fit_y = fitresult(down_fit_x);
down_fit_y = down_fit_y';

%fit up edge
[xData, yData] = prepareCurveData( up_errorX, up_errorY );
%ft = fittype( 'poly2' );
%[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 6.143039265005923E-7;
[fitresult, gof] = fit( xData, yData, ft, opts );
%find where fit reaches zero
I = up_max;
M = 1;
while M > mean_error
    M = fitresult(x(I));
    I = I-1;
end
up_fit_bot = I;
%make down edge fit line
up_fit_x = x(up_fit_bot:up_max);
up_fit_y = fitresult(up_fit_x);
up_fit_y = up_fit_y';

%avg height fit of flat part - moved up
%flat_error = y(down_max:up_min);
%mean_error = mean(flat_error);
%flat_fit_x = x(down_min:up_max);
%flat_fit_y = ones(length(flat_fit_x),1);
%flat_fit_y = flat_fit_y*mean_error;

%find intersections
x1 = flat_fit_x;
y1 = flat_fit_y;
x2 = down_fit_x;
y2 = down_fit_y;
[x0,y0] = intersections(x1,y1,x2,y2);
down_intersect = [x0,y0];
X1 = x0;
x2 = up_fit_x;
y2 = up_fit_y;
[x0,y0] = intersections(x1,y1,x2,y2);
up_intersect = [x0,y0];
X2 = x0;

%find diameter
diameter = X2-X1;
diameter_pixels = diameter/10;
diameter_odd_even = diameter_pixels;
edges = [X1 X2];

%figure;plot(x,y_deriv)
%hold on
%plot(x(down_min),y_deriv(down_min),'o',x(down_max),y_deriv(down_max),'o',x(up_min),y_deriv(up_min),'o',x(up_max),y_deriv(up_max),'o')
%hold off
tit1 = strcat('Error; file:',name,'; comp#:',comparison_number);
figure;plot(x,y,x,y_fit);title(tit1);
hold on
plot(down_fit_x,down_fit_y,up_fit_x,up_fit_y)
plot(flat_fit_x,flat_fit_y)
plot(down_intersect(1),down_intersect(2),'o',up_intersect(1),up_intersect(2),'o')
hold off
%IMG = Ipack.crop_sides;
%[a,b] = size(IMG);
%X1 = round(X1/10);
%X2 = round(X2/10);
%figure;imshow(IMG);line([1,b],[X1,X1],'color','r');line([1,b],[X2,X2],'color','r');

end