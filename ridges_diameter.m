function Ipack = ridges_diameter(Ipack)

polyDeg = Ipack.polynomialDegree;
slopeFilter = Ipack.slopeFilter;
flatFilter = Ipack.flatFilter;
ridgeArray = Ipack.crop_ridge;

[~,~,c] = size(ridgeArray);
ridge = cell(1,c);
for n = 1:c
    ridge{n} = ridgeArray(:,:,n);
end
down = cell(1,c);
up = cell(1,c);
flat = cell(1,c);
height = cell(1,c);
intersects = cell(1,c);
separations = cell(1,c);
diameterIntersection_pixels = cell(1,c);
diameterHeight = cell(1,c);

for n = 1:c
    ridge1 = ridge{n};
    [fitfringe, fitderiv] = estFit(ridge1);
    downridge = downRidge(ridge1,fitfringe, fitderiv, slopeFilter);
    upridge = upRidge(ridge1,fitfringe, fitderiv, slopeFilter);
    [flatridge,flatavg] = flatRidge(ridge1,fitfringe, fitderiv, flatFilter,upridge,downridge);
    [intersect,diameterIntersect_pixels] = diameterIntersections(ridge1,downridge,upridge,flatridge,polyDeg);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [diameterMaxHt,separation] = diameterMaxHeight(ridge1,flatridge,upridge,downridge);
    
    down{n}=downridge;
    up{n} = upridge;
    flat{n} = flatridge;
    height{n} = flatavg;
    intersects{n} = intersect;
    diameterIntersection_pixels{n} = diameterIntersect_pixels;
    separations{n}=separation;
    diameterHeight{n} = diameterMaxHt;
    
end

Ipack.downridge=down;
Ipack.upridge=up;
Ipack.flatridge=flat;
Ipack.fringeheight=height;
Ipack.intersects=intersects;
Ipack.separations=separations;
Ipack.diameterIntersection_pixels=diameterIntersection_pixels;
Ipack.diameterMaxHeight = diameterHeight;
Ipack.ridge=ridge;


end