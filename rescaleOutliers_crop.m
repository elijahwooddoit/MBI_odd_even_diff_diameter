function Ipack = rescaleOutliers_crop(Ipack)

%finds intensity outliers by comparing to average nearest neighbor pixel intesnsity
%and then replaces outlier intensity with nearest neighbor average if it is
%greater than 2*avgerage

int2 = Ipack.crop;
[p,q] = size(int2);
theChosen = [];

%find the values of the nearest neighbor points, checking first to make
%sure they exist, for every pixel in image int2
for w = 1:p
    for v = 1:q
        if (w+1) <= p
            up = int2((w+1),v);
        else
            up=nan;
        end
        %%%%
        if (w-1) >= 1
            down = int2((w-1),v);
        else
            down=nan;
        end
        %%%%
        if v+1 <= q
            right = int2(w,(v+1));
        else
            right=nan;
        end
        %%%%
        if v-1 >= 1
            left = int2(w,(v-1));
        else
            left = nan;
        end
        %%%%
        if (w+1) <= p && v+1 <= q
            upright = int2((w+1),(v+1));
        else
            upright = nan;
        end
        %%%%
        if (w+1) <= p && v-1 >= 1
            upleft = int2((w+1),(v-1));
        else
            upleft = nan;
        end
        %%%%
        if (w-1) >= 1 && v+1 <= q
            downright = int2((w-1),(v+1));
        else
            downright = nan;
        end
        %%%%
        if (w-1) >= 1 && v-1 >= 1
            downleft = int2((w-1),(v-1));
        else
            downleft = nan;
        end
        center = int2(w,v); %fill matrix of pixels to be replaced, named 'theChosen'
        neighbors = [up,upright,right,downright,down,downleft,left,upleft];
        neighbors(isnan(neighbors)) = [];
        avgI = mean(neighbors);
        if center > 2*avgI
            theChosen = [theChosen;w,v,avgI,center];
        elseif center > (avgI+0.15)
            theChosen = [theChosen;w,v,avgI,center];
        elseif center < (avgI-0.15)
            theChosen = [theChosen;w,v,avgI,center];
        end
        %%%%
    end
end

%replace stray pixels defined in theChosen
int3 = int2;
[s,t] = size(theChosen);
for w = 1:s
    int3(theChosen(w,1),theChosen(w,2))=theChosen(w,3);
end

%rescale image int3 and save to Ipack
int3 = int3-min(min(int3));
x = 1 / max(max(int3));
int3 = int3 * x;
Ipack.crop_rescaled = int3;
end