function [flatridge,flatavg] = flatRidge(ridge1,fitfringe, fitderiv, flatFilter,upridge,downridge)

x=ridge1(:,1);
y=ridge1(:,2);

estFit = fitfringe;
estFitDeriv = fitderiv;
x1 = estFit(:,1);
y1 = estFit(:,2);
x2 = estFitDeriv(:,1);
y2 = estFitDeriv(:,2);

X2 = x2;
Y2 = y2;
n = 1;
while n <= length(Y2)   %finds all values with slope greater than / equal to filter defined in master file and puts them in up
    if Y2(n) > flatFilter || Y2(n) < (flatFilter*(-1))
        Y2(n) = [];
        X2(n) = [];
    else
        n=n+1;
    end
end
flat = [X2,Y2];
flatx = flat(:,1);

minVal = downridge(end,1);
maxVal = upridge(1,1);

n=1;
while n <= length(flatx)   %finds any flat curve outside beyond the fringe edges and deletes it
    if flatx(n) > maxVal || flatx(n) < minVal
        flatx(n) = [];
        flat(n,:) = [];
    else
        n=n+1;
    end
end

flatridge = [];
for n = 1:length(x) %uses down as a template to find the actual fringe values spanning the same range
    if x(n) >= flat(1,1) && x(n) <= flat(end,1)
        flatridge = [flatridge; x(n), y(n)];
    end
end

figure; plot(flatridge(:,1),flatridge(:,2))
%Ipack.flatridge = flatridge;

figure; plot(flatridge(:,1),flatridge(:,2),upridge(:,1),upridge(:,2),downridge(:,1),downridge(:,2))

flatavg = mean(flatridge(:,2)); %calculate average height of fringe so that they can be lined up in future
%Ipack.flatavg = flatavg;

end