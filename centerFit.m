function [cent,deltac] = centerFit(static,dynam,initbounds)

%static
%dynam
%initbounds = fitbounds{2}; %%%%%

[j,k] = size(dynam);
[h,i] = size(static);
l = initbounds(1);
u = initbounds(2);
c = initbounds(3);
dynam1 = dynam(:,(c-l):(c+u));


err = zeros(h,i);
cspan = (l+1):(k-u-1);
error1 = zeros(length(cspan),2);
for n = 1:length(cspan)
    c1 = cspan(n);
    dynam1 = dynam(:,(c1-l):(c1+u));
    for m = 1:h
        for w = 1:i
            err(m,w) = (static(m,w)-dynam1(m,w))^2;
        end
    end
    error = sum(sum(err));
    error1(n,1) = c1;
    error1(n,2) = error;
end

figure;plot(error1(:,1),error1(:,2))
[f,g] = min(error1(:,2));
cent = error1(g,1);
deltac = cent-c;

end