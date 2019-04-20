function upridge = upRidge(ridge1,fitfringe, fitderiv, slopeFilter)

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
    if Y2(n) < slopeFilter
        Y2(n) = [];
        X2(n) = [];
    else
        n=n+1;
    end
end
up = [X2,Y2];

upx = up(:,1);
A=[];
B=[];
for n = 2:length(up)  %finds all separations in n greater than 1 to find discontinuous segments
    if (upx(n)-upx(n-1)) ~= 1
        A = [A; upx(n)];
        B = [B; n];
    end
end

C=[];
if isempty(A)   %if discontinuous segments were found decides which is the right one by finding the biggest segment
else
    upx=[upx;0];
    B = [1;B;length(upx)];
    for n = 2:length(B) %finds the size of the segments and puts them in vector C
            temp = upx(B(n)-1)-upx(B(n-1));
            C = [C; temp];
    end
    upx(end)=[];
    C = C+1;
    B(end)=[length(upx)];
    
    [M,I] = max(C); %finds the largest segment, makes that up and deletes the rest
    I = I+1;
    temp = [B(I-1),(B(I)-1)];
    up = up(temp(1):temp(2),:);
end

if length(up) > 50    %if down has more than 50 elements it is shortened to the inner 50 elements
    up = up(1:50,:);
end

upridge = [];
for n = 1:length(x) %uses down as a template to find the actual fringe values spanning the same range
    if x(n) >= up(1,1) && x(n) <= up(end,1)
        upridge = [upridge; x(n), y(n)];
    end
end

figure; plot(upridge(:,1),upridge(:,2))
%Ipack.upridge = upridge;

end