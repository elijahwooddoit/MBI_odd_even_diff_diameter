function Ipack = crop7_Ipack(Ipack)


%peakProminence = 0.2;
%peakDist = 10;

I = Ipack.crop;
I2 = sum(I,1);

I2 = I2-min(min(I2));       %I2 is intensity integration of cropped file
o = 1 / max(max(I2));
I2 = I2 * o;

[~,locs]=findpeaks(I2,'MinPeakProminence',0.2,'MinPeakDistance',10);  %find intensity peaks in I2 spectrum

fringeCenters = zeros(length(locs),3);   %define center area as +/- 40 from avg center
fringeCenters(:,2) = locs;
fringeCenters(:,1) = locs-40;
fringeCenters(:,3) = locs+40;


%%% important! need to figure out how to change the mumbo-jumbo below to
%%% something that makes a slice of the fringe image
Ifringe = cell(1,length(locs));
center_intensities = cell(1,length(locs));
crops = zeros(length(locs),2);
for n = 1:length(locs)
    %integrate the total intensities over the erea of where the fringes
    %should be and filter the resulting intensity vector
    Ifringe{n} = I(:,fringeCenters(n,1):fringeCenters(n,3));
    
    center_intensities{n} = sum(Ifringe{n},2);
    y = sgolayfilt(center_intensities{n},2,11);
    %smooth out the intensity  plot
    y4 = sgolayfilt(center_intensities{n},1,41);
    y4 = smooth(y4,'lowess');
    x = 1:length(y);
    %figure;plot(x,y)
    [xData, yData] = prepareCurveData( x, y4 );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    y_compare = fitresult(x);
    y4 = y-y_compare;
    y5 = gradient(y4);
    %threshold the intensity vector to 1's and zeros
    y2 = 1:length(y);
    threshold = mean(y4);
    y2(find(y4 < threshold)) = 0;
    y2(find(y4 >= threshold)) = 1;
    %split the intensity vector into several vectors separated by zeros and
    %find the largest one returned in l
    v = y2;
    w = [false v~=0 false]; %// "close" v with zeros, and transform to logical
    starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
    ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
    result = arrayfun(@(s,e) v(s:e), starts, ends, 'uniformout', false); %// build result
    [~,d] = size(result);
    l = 1;
    k = 1;
    for m = 1:d
        if length(result{m}) > k
            k = length(result{m});
            l = m;
        end
    end
    %find the indices of the largest peak
    %plot(y4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = [0,y2];
    [~,positions] = findpeaks(y3);
    peakstart = positions(l) - 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Index exceeds matrix dimensions.
    o = 1;
    p = 0;
    while o == 1
        p = p+1;
        o = y2(peakstart+p);
    end
    peakend = peakstart + p - 1;
    %smooth out the intensity  plot - moved up to beginning of loop
    %y4 = sgolayfilt(center_intensities{n},1,41);
    %y4 = smooth(y4,'lowess');
    %y5 = gradient(y4);
    
    %find the upper and lower limits to crop the image vertically
    o = 1;
    p = 0;
    while o >= 0
        p = p+1;
        o = y5(peakstart-p);
        if (peakstart-p) == 1
            o = -1;
        end
    end
    upcrop = peakstart-p;
    o = -1;
    p = 0;
    while o <= 0
        p = p+1;
        o = y5(peakend+p);
        if (peakend+p) == length(y5)
            o = -1;
        end
    end
    downcrop = peakend+p;
    crops(n,1) = upcrop;
    crops(n,2) = downcrop;
end

upcrop = mean(crops(:,1));
upcrop = round(upcrop);
downcrop = mean(crops(:,2));
downcrop = round(downcrop);
Icrop = I(upcrop:downcrop,:);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Warning: Integer operands are required for colon operator when used as index
Ipack.crop = Icrop;


%visualization
%x = 1:length(I2);
%figure; plot(x,I2,locs,pks,'o')
%figure;plot(y)
%figure;plot(y2)
%figure;plot(y4)
%figure;plot(y5)
%figure;imshow(I)
figure;imshow(Icrop)

end