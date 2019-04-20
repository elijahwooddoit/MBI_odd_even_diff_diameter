function fringeCenter = fitFringeCenter(center, I, I2, Ipack)

window = 20;
peakProminence = Ipack.peakProminence;
peakDist = Ipack.peakDist;

if (center - window) > 0    %define area of interest as the intensity peak plus a window to either side, check to make sure indices are ok
    lowerbound = center - window;
else
    lowerbound = 1;
end
if (center + window) <= length(I2)
    upperbound = center + window;
else
    upperbound = length(I2);
end
I3 = I(:,lowerbound:upperbound);

dummyFringe = zeros(1,size(I3,1)); %%%find the number of peaks within the vertical limits and interpret this as the number of fringes
dummyFringe = dummyFringe';
for w=1:length(dummyFringe)
    A = I3(w,:);
    [pks,locs1]=findpeaks(A,'MinPeakProminence',peakProminence,'MinPeakDistance',peakDist);
    if isempty(locs1)
    else
        [~,k] = max(pks);
        k = k(1);
        dummyFringe(w) = locs1(k);
    end
end
x1 = 1:length(dummyFringe);
x1 = x1';
dummyFringe1 = [x1, dummyFringe];

[c, ~] = find(dummyFringe == 0); %finds any rows that still have dummy zero values and eliminates them
dummyFringe1(c,:) = [];


%figure;plot(dummyFringe1(:,1),dummyFringe1(:,2))

x2 = dummyFringe1(:,1); % fits data within window w/ 4th order polynomial and finds the derivative of the fit
y2 = dummyFringe1(:,2);
[xData, yData] = prepareCurveData( x2, y2 );
ft = fittype( 'poly4' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'LAR';
[fitresult, gof] = fit( xData, yData, ft, opts );
x3 = min(x2):max(x2);
x3 = x3';
y = fitresult(x3);
x4 = x3;
x4(length(x3)) = [];
y1 = diff(y);
%figure;plot(x2,y2,'o',x3,y)
%figure;plot(x4,y1)

flatFilter = 0.12;
X2 = x4;
Y2 = y1;
n = 1;
while n <= length(Y2)   %finds all values with slope greater than / equal to flat filter deletes them
    if Y2(n) > flatFilter || Y2(n) < (flatFilter*(-1))
        Y2(n) = [];
        X2(n) = [];
    else
        n=n+1;
    end
end
flat = [X2,Y2];
flatx = flat(:,1);
%plot(flat(:,1),flat(:,2))
fringeCenter = mean(X2);


end