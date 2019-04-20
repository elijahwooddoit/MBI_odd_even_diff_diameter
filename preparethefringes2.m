function Ipack = preparethefringes2(Ipack)

flatridge = Ipack.flatridge2;
upridge = Ipack.upridge2;
downridge = Ipack.downridge2;
BoundsBase = Ipack.boundsBase2;
BoundsFit = Ipack.boundsFit2;
lowhigh = Ipack.fitlowhigh2;
lowhighBase = Ipack.baselowhigh2;
Ifilt = Ipack.ImEnhance;
[~,b] = size(flatridge);
lb = lowhighBase(1);
hb = lowhighBase(2);
lf = lowhigh(1);
hf = lowhigh(2);

fitbounds = cell(1,b);
Iflat = cell(1,b);  %creates cropped fringe images Iflat from the larger image and defines bounds within this image
for n = 1:b
    Iflat{n} = Ifilt(BoundsBase(n,4):BoundsBase(n,5),BoundsBase(n,1):BoundsBase(n,2));
    fitbounds{n} = [lowhigh(1),lowhigh(2),(BoundsBase(n,3)-BoundsBase(n,1)+1),BoundsBase(n,4),BoundsBase(n,5)];%%for fit in base Iflat: dist left, dist right, center, top, bottom
end

Iregr1 = Ifilt(BoundsBase(1,4):BoundsBase(1,5),(BoundsBase(1,3)-lowhigh(1)):(BoundsBase(1,3)+lowhigh(2)));
Bounds = BoundsBase;    %%%%or BoundsFit
for n = 2:b
    static = Iregr1;
    dynam = Iflat{n};
    initbounds = fitbounds{n};
    [cent,deltac] = centerFit(static,dynam,initbounds);
    Bounds(n,3) = Bounds(n,3)+deltac;
    Bounds(n,1) = Bounds(n,3)-lb;   %%%%or lf
    Bounds(n,2) = Bounds(n,3)+hb;   %%%%or hf
end

%need to do another loop to figure out top and bottom - start of edges

[~,u] = size(upridge);
uplim = zeros(u,1);
for n = 1:u
    uplim(n) = min(upridge{n}(:,1));
end

[~,u] = size(downridge);
downlim = zeros(u,1);
for n = 1:u
    downlim(n) = max(downridge{n}(:,1));
end

Bounds(:,4) = downlim;
Bounds(:,5) = uplim;
temp1=min(Bounds(:,4));
temp2=max(Bounds(:,5));
temp11=(ones(u,1))*temp1;
temp22=(ones(u,1))*temp2;
Bounds(:,4) = temp11;
Bounds(:,5) = temp22;

sizecheck = zeros(u,2);
for n = 1:u
    sizecheck(n,1)=Bounds(n,2)-Bounds(n,1)+1;
    sizecheck(n,2)=Bounds(n,5)-Bounds(n,4)+1;
end

Ipack.bounds = Bounds; %non-cell matrix, column 1 = left boundar, column 2 = right boundary, column 3 = fringe center, column 4 = lower boundary, column 5 = upper boundary
Ipack.szcheck = sizecheck;


end