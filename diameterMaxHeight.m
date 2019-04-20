function [diameterMaxHt,separation] = diameterMaxHeight(ridge1,flatridge,upridge,downridge)

heightFilter = 1.5;

lowerbound = mean(downridge(:,1));
upperbound = mean(upridge(:,1));
ridge2 = ridge1;
ridge2x = ridge2(:,1);
n = 1;
while n <= length(ridge2)   %only keeps ridge values up to halfway up fringe edges
    if ridge2(n) < lowerbound || ridge2(n) > upperbound
        ridge2(n,:) = [];
    else
        n=n+1;
    end
end


lowerbound = min(downridge(:,1));   %fit flat part of fringe to line
upperbound = max(upridge(:,1));
flatx = flatridge(:,1);
flaty = flatridge(:,2);
[xData, yData] = prepareCurveData( flatx, flaty );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
x1 = lowerbound:upperbound;
x1 = x1';
y1 = fitresult(x1);
y2 = y1+heightFilter;

%plot(ridge2(:,1),ridge2(:,2),x1,y1)

[xData, yData] = prepareCurveData( ridge2(:,1), ridge2(:,2) );  %fit smoothing curve to fringe with edges
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
opts.SmoothingParam = 0.9999;
[fitresult, gof] = fit( xData, yData, ft, opts );
x = min(ridge2(:,1)):max(ridge2(:,1));
x=x';
y = fitresult(x);

[x0,y0] = intersections(x1,y2,x,y);
[~,minx] = min(x0);
[~,maxx] = max(x0);
separation = [x0(minx),y0(minx);x0(maxx),y0(maxx)];

plot(ridge2(:,1),ridge2(:,2),'o',x1,y1,x,y,separation(:,1),separation(:,2),'*')

diameterMaxHt = separation(2,1)-separation(1,1);

end

