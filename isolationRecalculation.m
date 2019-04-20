function Ipack = isolationRecalculation(Ipack)

BoundsBase = Ipack.boundsBase;
BoundsFit = Ipack.boundsFit;
lowhigh = Ipack.fitlowhigh;
lowhighBase = Ipack.baselowhigh;
upridge = Ipack.upridge;
downridge = Ipack.downridge;
flatridge = Ipack.flatridge;
ridge = Ipack.ridge;
Ifilt = Ipack.crop_smoothed2;
Ienhance = Ipack.ImEnhance;

%Ibase = Ifilt(BoundsBase(1,4):BoundsBase(1,5),BoundsBase(1,1):BoundsBase(1,2));
%Ifit = Ifilt(BoundsFit(1,4):BoundsFit(1,5),BoundsFit(1,1):BoundsFit(1,2));
BoundsBase2 = (BoundsBase-1).*10+1;
BoundsFit2 = (BoundsFit-1).*10+1;
I4 = Ienhance(BoundsFit2(1,4):BoundsFit2(1,5),BoundsFit2(1,1):BoundsFit2(1,2));

lowhigh2 = lowhigh.*10;
lowhighBase2 = lowhighBase.*10;

upridge2 = upridge;
[~,b] = size(upridge);
for n=1:b
    upridge2{n} = (upridge{n}-1).*10+1;
end

downridge2 = downridge;
[~,b] = size(downridge);
for n=1:b
    downridge2{n} = (downridge{n}-1).*10+1;
end

flatridge2 = flatridge;
[~,b] = size(flatridge);
for n=1:b
    flatridge2{n} = (flatridge{n}-1).*10+1;
end

ridge2 = ridge;
[~,b] = size(ridge);
for n=1:b
    ridge2{n} = (ridge{n}-1).*10+1;
end

%z = zeros(length(flatridge2{1}),1);
%for n=1:length(flatridge2{1})
%    z(n) = Ienhance(flatridge2{1}(n,1),flatridge2{1}(n,2));
%end
%ridgeline = [flatridge2{1},z];
%[YY,XX] = size(Ienhance);
%yy = 1:YY;
%xx = 1:XX;
%[xxx,yyy]=meshgrid(xx,yy);
%figure;surf(xxx,yyy,Ienhance,'LineStyle','none');line(ridgeline(:,2),ridgeline(:,1),ridgeline(:,3),'color','r')

Ipack.boundsBase2 = BoundsBase2;
Ipack.boundsFit2 = BoundsFit2;
Ipack.fitlowhigh2 = lowhigh2;
Ipack.baselowhigh2 = lowhighBase2;
Ipack.upridge2 = upridge2;
Ipack.downridge2 = downridge2;
Ipack.flatridge2 = flatridge2;
Ipack.ridge2 = ridge2;


end