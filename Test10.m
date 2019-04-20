error_data = error_plots{2};
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
%
figure;plot(x,y,x,y_fit)
figure;plot(x,y_deriv)
