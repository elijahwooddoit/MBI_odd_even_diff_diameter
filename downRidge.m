function downridge = downRidge(ridge1,fitfringe, fitderiv, slopeFilter)

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
while n <= length(Y2)   %finds all values with slope less than / equal to filter defined in master file and puts them in down
    if Y2(n) > (slopeFilter*(-1))
        Y2(n) = [];
        X2(n) = [];
    else
        n=n+1;
    end
end
down = [X2,Y2];

%downx = down(:,1);
%A=[];
%B=[];
%for n = 2:length(down)  %finds all separations in n greater than 1 to find discontinuous segments
%    if (downx(n)-downx(n-1)) ~= 1
%        A = [A; downx(n)];
%        B = [B; n];
%    end
%end

%%%%%%%%%%%%
%if isempty(A)   %if a discontinuous segment has been found, only the first segment is kept
%else
%    n = A(1);
%    w = find(downx == n);
%    down(w(1):end,:)=[];
%end

downx = down(:,1);
A=[];
B=[];
for n = 2:length(down)  %finds all separations in n greater than 1 to find discontinuous segments
    if (downx(n)-downx(n-1)) ~= 1
        A = [A; downx(n)];
        B = [B; n];
    end
end

C=[];
if isempty(A)   %if discontinuous segments were found decides which is the right one by finding the biggest segment
else
    downx=[downx;0];
    B = [1;B;length(downx)];
    for n = 2:length(B) %finds the size of the segments and puts them in vector C
            temp = downx(B(n)-1)-downx(B(n-1));
            C = [C; temp];
    end
    downx(end)=[];
    C = C+1;
    B(end)=[length(downx)];
    
    [M,I] = max(C); %finds the largest segment, makes that down and deletes the rest
    I = I+1;
    temp = [B(I-1),(B(I)-1)];
    down = down(temp(1):temp(2),:);
end

if length(down) > 50    %if down has more than 50 elements it is shortened to the inner 50 elements
    down = down(end-49:end,:);
end

downridge = [];
for n = 1:length(x) %uses down as a template to find the actual fringe values spanning the same range
    if x(n) >= down(1,1) && x(n) <= down(end,1)
        downridge = [downridge; x(n), y(n)];
    end
end

figure; plot(downridge(:,1),downridge(:,2))
%Ipack.downridge = downridge;

end