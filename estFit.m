function [fitfringe, fitderiv] = estFit(ridge1)

x=ridge1(:,1);
y=ridge1(:,2);

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'poly6' );

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult, xData, yData );

x1 = min(x):max(x);
x1 = x1';
y1 = fitresult(x1);

x2 = x1;
x2(length(x1)) = [];
y2 = diff(y1);

%x3 = x2;
%x3(length(x2)) = [];
%y3 = diff(y2);

figure;plot(x,y,'o',x1,y1)
figure;plot(x2,y2)
%figure;plot(x3,y3)

fitfringe = [x1,y1];
fitderiv = [x2,y2];

end