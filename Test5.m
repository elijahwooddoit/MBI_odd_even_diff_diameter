I = Ipack.crop_smoothed2;
[a,b]=size(I);
error = zeros(a,1);
for n = 1:a
    t = gradient(I(n,:));
    tt = sgolayfilt(t,2,41);
    w = gradient(tt);
    ww = sgolayfilt(w,2,41);
    w_inv = ww * -1;    %max curvature ridge line
    r = abs(tt);
    rr = sgolayfilt(r,2,41);    %ridge edges max sloping curvature
    rr = rr - mean(rr);
    k = find(rr < 0);
    rr(k) = 0;
    rrr = sgolayfilt(rr,2,11);
    threshold = max(rrr)/10;
    [~,locs1] = findpeaks(rrr,'MinPeakProminence',threshold);
    q = -1 * rr;
    q = q+(0-min(q));
    threshold = min(ww)*0.3;
    k = find(ww < threshold);
    q(k) = q(k)*10;
    qq = sgolayfilt(q,2,41);    %zero slope ridge line
    qq = qq - mean(qq);
    k = find(qq < 0);
    qq(k) = 0;
    qqq = sgolayfilt(qq,2,21);
    threshold = max(qqq)/10;
    [~,locs2] = findpeaks(qqq,'MinPeakProminence',threshold);
    locs = locs2;
    N = 84;
    dist = abs(locs-N);
    minerror = min(dist);
    error(n) = minerror;
end

%error = error-min(error);
%temp = 1/max(error);
%error = error * temp;
figure;plot(error)
error1 = sgolayfilt(error,2,51);
figure;plot(error1)




