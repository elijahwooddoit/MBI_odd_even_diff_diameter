function Ipack = isolate3Dfringe(Ipack)
%%
height = Ipack.fringeheight;
Im = Ipack.crop_smoothed;

flatridge = Ipack.flatridge;    %repeat the process of finding the shortest fringe and shortening the other fringes
[~,b] = size(flatridge);
dummy = zeros(b,1);
for n = 1:b
    a = length(flatridge{n});
    dummy(n) = a;
end
[M,I] = min(dummy);
span1 = [min(flatridge{I}(:,1)),max(flatridge{I}(:,1))];
lg = length(flatridge{I});
for n = 1:b
    [C,ia,ib] = intersect(flatridge{I}(:,1),flatridge{n}(:,1));
    ib = ib';
    flatridge{n} = flatridge{n}(ib,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ln1 = flatridge{1};
ln2 = flatridge{2};
[p,q]=size(Im);
P=1:p;
Q=1:q;
[pp,qq]=meshgrid(Q,P);
figure;surf(pp,qq,Im);line(ln1(:,2),ln1(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

borders = cell(1,b);
for n = 1:b
    borders{n} = zeros(length(flatridge{n}),4); %borders first column is row in image, second column is distance to front, third column is ridge line, fourth column is distannce to back column in image
    for m=1:length(flatridge{n})    %creating smaller image of just fringe
        u = flatridge{n}(m,1);
        borders{n}(m,1)=u;
        w = flatridge{n}(m,2);
        w1 = w;
        borders{n}(m,3)=w1;
        v = Im(u,w1);
        v1 = v;
        t = 0;
        while v > 0.2*v1
            try
               w = w+1;
               v = Im(u,w);
               t = t+1;
            catch
                disp(v);disp(w);disp(t);
            end
            %w = w+1;
            %v = Im(u,w);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %t = t+1;
        end
        borders{n}(m,4) = w;
        w = w1;
        v = v1;
        t = 0;
        while v > 0.2*v1
            w = w-1;
            v = Im(u,w);
            t = t+1;
        end
        borders{n}(m,2) = w;
    end
end

MinMidMax = zeros(b,5);
for n = 1:b
    limits = borders{n};
    minlim = min(limits(:,2));
    maxlim = max(limits(:,4));
    fringe = round(mean(limits(:,3)));
    MinMidMax(n,1) = minlim;
    MinMidMax(n,3) = maxlim;
    MinMidMax(n,2) = fringe;
    MinMidMax(n,4) = fringe - minlim;
    MinMidMax(n,5) = maxlim - fringe;
end

low = max(MinMidMax(:,4));
high = max(MinMidMax(:,5));
span2 = low+high+1;
lowhighBase = [low, high];

minmax = zeros(b,5);
for n = 1:b
    minmax(n,1) = MinMidMax(n,2)-low;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minmax(n,2) = MinMidMax(n,2)+high;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minmax(n,3) = MinMidMax(n,2);
    minmax(n,4) = span1(1);
    minmax(n,5) = span1(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write a loop that prevents cropping indices (minmax columns 1 and 2) lower
%than one or above the image size and fills out what would be the rest of
%the image with whatever the minimum value in that cropping is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Im2 = cell(1,b);
for n = 1:b
    Im2{n}=Im(span1(1):span1(2),minmax(n,1):minmax(n,2));
end

f = span1(1):span1(2);
span = length(f);
X=1:span2;
Y=1:span;
[xx,yy]=meshgrid(X,Y);
%for n = 1:b
%    J = Im2{n};
%    figure;surf(xx,yy,J);
%end

Ibase = Im2;    %cell that contains base images of the fringes
BoundsBase = minmax;    %box boundaries within image corresponding to Ibase

%%
height = Ipack.fringeheight;
Im = Ipack.crop_smoothed;

flatridge = Ipack.flatridge;    %repeat the process of finding the shortest fringe and shortening the other fringes
[~,b] = size(flatridge);
dummy = zeros(b,1);
for n = 1:b
    a = length(flatridge{n});
    dummy(n) = a;
end
[M,I] = min(dummy);
span1 = [min(flatridge{I}(:,1)),max(flatridge{I}(:,1))];
lg = length(flatridge{I});
for n = 1:b
    [C,ia,ib] = intersect(flatridge{I}(:,1),flatridge{n}(:,1));
    ib = ib';
    flatridge{n} = flatridge{n}(ib,:);
end

borders = cell(1,b);
for n = 1:b
    borders{n} = zeros(length(flatridge{n}),4); %borders first column is row in image, second column is distance to front, third column is ridge line, fourth column is distannce to back column in image
    for m=1:length(flatridge{n})    %creating smaller image of just fringe
        u = flatridge{n}(m,1);
        borders{n}(m,1)=u;
        w = flatridge{n}(m,2);
        w1 = w;
        borders{n}(m,3)=w1;
        v = Im(u,w1);
        v1 = v;
        t = 0;
        while v > 0.4*v1
            w = w+1;
            v = Im(u,w);
            t = t+1;
        end
        borders{n}(m,4) = w;
        w = w1;
        v = v1;
        t = 0;
        while v > 0.4*v1
            w = w-1;
            v = Im(u,w);
            t = t+1;
        end
        borders{n}(m,2) = w;
    end
end

MinMidMax = zeros(b,5);
for n = 1:b
    limits = borders{n};
    minlim = min(limits(:,2));
    maxlim = max(limits(:,4));
    fringe = round(mean(limits(:,3)));
    MinMidMax(n,1) = minlim;
    MinMidMax(n,3) = maxlim;
    MinMidMax(n,2) = fringe;
    MinMidMax(n,4) = fringe - minlim;
    MinMidMax(n,5) = maxlim - fringe;
end

low = max(MinMidMax(:,4));
high = max(MinMidMax(:,5));
span2 = low+high+1;

minmax = zeros(b,5);
for n = 1:b
    minmax(n,1) = MinMidMax(n,2)-low;
    minmax(n,2) = MinMidMax(n,2)+high;
    minmax(n,3) = MinMidMax(n,2);
    minmax(n,4) = span1(1);
    minmax(n,5) = span1(2);
end

Im2 = cell(1,b);
for n = 1:b
    Im2{n}=Im(span1(1):span1(2),minmax(n,1):minmax(n,2));
end

f = span1(1):span1(2);
span = length(f);
X=1:span2;
Y=1:span;
[xx,yy]=meshgrid(X,Y);
%for n = 1:b
%    J = Im2{n};
%    figure;surf(xx,yy,J);
%end

Ifit = Im2;    %cell that contains base images of the fringes
BoundsFit = minmax;    %box boundaries within image corresponding to Ibase
lowhigh = [low,high];

Ipack.Ibase = Ibase;
Ipack.boundsBase = BoundsBase;
Ipack.boundsFit = BoundsFit;
Ipack.fitlowhigh = lowhigh;
Ipack.baselowhigh = lowhighBase;
end