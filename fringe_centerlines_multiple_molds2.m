function Ipack = fringe_centerlines_multiple_molds2(Ipack)

upridge = Ipack.upridge2;
downridge = Ipack.downridge2;
flatridge = Ipack.flatridge2;
BoundsBase = Ipack.boundsBase2;
BoundsFit = Ipack.boundsFit2;
bounds = Ipack.bounds;
lowhigh = Ipack.fitlowhigh2;
lowhighBase = Ipack.baselowhigh2;
ridge = Ipack.ridge2;
I = Ipack.ImEnhance;
[~,b] = size(flatridge);

%%%% shorten ridge using up and downridge (inner 25 of each) before it goes
%%%% into the border finding loop
%downtemp = downridge{1}(1:end-25,1);
%uptemp = upridge{1}(1:25,1);
%for m = 1:b
ridgemin = min(downridge{1}(:,1));
ridgemax = max(upridge{1}(:,1));
for m = 1:b    
    n = 1;
    while n <= length(ridge{m})
        if ridge{m}(n,1) < ridgemin
            ridge{m}(n,:) = [];
        elseif ridge{m}(n,1) > ridgemax
            ridge{m}(n,:) = [];
        else
            n = n+1;
        end
    end
end

% create an interpolation fit of the fringe with only integer points to use
% in finding the fringe edges
x1 = cell(1,b);
y1 = cell(1,b);
for n = 1:b
    x = ridge{n}(:,1);
    y = ridge{n}(:,2);
    [xData, yData] = prepareCurveData( x, y );
    ft = 'pchipinterp';
    [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    x1{n} = min(x):max(x);
    y1{n} = fitresult(x1{n});
    y1{n} = round(y1{n});
end


%borders = cell(1,b);    %column 1=pixel row, column 2 = fringe front, column 3 = ridgeline, column 4 = fringe back
%for n = 1:b
borders = cell(1,b);
filter = 0.8;
for n=1:b
    borders{n} = zeros(length(x1{n}),4);
    for m=1:length(x1{n})
        u = x1{n}(m);
        borders{n}(m,1)=u;
        w = y1{n}(m);
        w1 = w;
        borders{n}(m,3)=w1;
        v = I(u,w1);
        v1 = v;
        t = 0;
        while v > filter*v1
            w = w+1;
            v = I(u,w);
            t = t+1;
        end
        borders{n}(m,4) = w;
        w = w1;
        v = v1;
        t = 0;
        while v > filter*v1
            w = w-1;
            v = I(u,w);
            t = t+1;
        end
        borders{n}(m,2) = w;
    end
end

%%%% define borders based on the centerline in bounds
centers = bounds(:,3);
dim1 = cell(1,b);
dim2 = cell(1,b);
for m = 1:b
    dim1{m} = borders{m}(:,1);
    dim2{m} = borders{m}(:,2:4);
    for n = 1:length(dim1{m})
        dim2{m}(n,:) = dim2{m}(n,:)-centers(m);
    end
end

%find largest dimensions of all the fringes to define box sizes
min_dim = zeros(1,b);
max_dim = zeros(1,b);
for n =1:b
    min_dim(n) = min(min(dim2{n}));
    max_dim(n) = max(max(dim2{n}))-min(min(dim2{n}))+1;
end

boxsz1 = length(dim1{1});
boxsz2 = max(max_dim)-min(min_dim)+1;
minmin = min(min_dim);
Im = zeros(boxsz1,boxsz2);
fringeIm = cell(1,b);

%%%% create image files for each fringe with the borders filled in and the
%%%% rest of the matrix set to zero
for m = 1:b
    cent = centers(m);
    fringeIm{m} = Im;
    for n = 1:length(dim1{m})
        dim1n = dim1{m}(n);
        leftbound = dim2{m}(n,1);
        rightbound = dim2{m}(n,3);
        leftbound1 = cent+leftbound;    %positions in I
        rightbound1 = cent+rightbound;
        %pos_in_Im = x - minmin + 1
        lngth = rightbound-leftbound;
        pix_row = I(dim1n,leftbound1:rightbound1);
        first_pos = leftbound - minmin +1;
        last_pos = first_pos+lngth;
        fringeIm{m}(n,first_pos:last_pos) = pix_row;
    end
end

%%%% gonna try to make a wieghted average center
wieghted_center = cell(1,b);
wieghted_center_round = cell(1,b);
center_path = zeros(boxsz1,2);
center_path2 = zeros(boxsz1,3);
for n = 1:b
    wieghted_center{n} = center_path;
    wieghted_center_round{n} = center_path2;
    for m = 1:boxsz1
        pix_row = fringeIm{n}(m,:);
        mean_matrix = 1:length(pix_row);
        weighted_center{n}(m,2) = sum(mean_matrix.*pix_row) / sum(pix_row);
        weighted_center{n}(m,1) = dim1{n}(m);
        wieghted_center_round{n}(m,1) = m;
    end
    wieghted_center_round{n}(:,2) = round(weighted_center{n}(:,2));
    
    for m = 1:boxsz1
        wieghted_center_round{n}(m,3) = fringeIm{n}(m,wieghted_center_round{n}(m,2));
    end
end


%%%% also change shift the tops of the fringes all up to 1
fringeImShift = fringeIm;
for n = 1:b
    for m=1:boxsz1
        pix_row = fringeIm{n}(m,:);
        maxpt = max(pix_row);
        shift = 1-maxpt;
        %pix_row = pix_row+shift;
        for w = 1:length(pix_row)
            if pix_row(w) ~= 0
                pix_row(w) = pix_row(w)+shift;
            end
        end
        fringeImShift{n}(m,:) = pix_row;
    end
end

%%%% then just hijack the error calculation from test1
error2 = zeros(boxsz1,1);
for n = 2:b
    for m=1:boxsz1
        pix_row = fringeImShift{n}(m,:);
        pix_row1 = fringeImShift{1}(m,:);
        avgElements = (nnz(pix_row)+nnz(pix_row1))/2;
        error = 0;
        for w = 1:length(pix_row)
            error = error+(pix_row1(w)-pix_row(w))^2;
        end
        error2(m) = error/avgElements;
    end
end
fringe_error = error2;
temp = ridgemin:ridgemax;
temp = temp';
fringe_error = [temp, fringe_error];
Ipack.fringe_error = fringe_error;

%visualization
figure;plot(error2)
n=1;
I4 = fringeIm{n};
[YY,XX] = size(I4);
yy = 1:YY;
xx = 1:XX;
[xxx,yyy]=meshgrid(xx,yy);
figure;surf(xxx,yyy,I4,'LineStyle','none');line(wieghted_center_round{n}(:,2),wieghted_center_round{n}(:,1),wieghted_center_round{n}(:,3),'color','red')
%hold on
%n=2;
%I4 = fringeIm{n};
%line(wieghted_center_round{n}(:,2),wieghted_center_round{n}(:,1),wieghted_center_round{n}(:,3),'color','red')
%surf(xxx,yyy,I4,'LineStyle','none')
%hold off

end