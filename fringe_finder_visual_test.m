
%this script will plot the ridges over the 3D image

%open up the image and define its size
I = Ipack.crop_smoothed2;
[a,b] = size(I);
[x,y]=size(I);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
figure;surf(xx,yy,I);
%figure;surf(xx,yy,Ibounds2);

%create a cell with fringe arrays of each row and a vector for each row
%consisting of all the important fringe points and one for each of the
%important things
fringes = cell(1,a);
fringes2 = cell(1,a);
ridges = cell(1,a);
fringe_middles = cell(1,a);
outer_bounds = cell(1,a);
inner_bounds = cell(1,a);
for n = 1:a
    %disp(n)
    %store fringe array in fringes, which is just the collection of fringe
    %array matrices from every row
    fringe_array = fringe_finder(I,n);
    fringes{n} = fringe_array;
    %store only data about fringe locations (ridges, edfes, etc.) in
    %fringes2
    fringe_col = [1 2 4 5 6 7 8 9 10 11 12 13];
    fringe_array2 = fringe_array(:,fringe_col);
    fringe_array2 = fringe_array2';
    fringe_array2 = fringe_array2(:)';
    fringe_array2 = sort(fringe_array2);
    fringe_array2 = fringe_array2(fringe_array2 ~= 0);
    fringes2{n} = fringe_array2';
    %store only data about ridges from each row in the cell ridges
    fringe_col = [9 10 12 13];
    ridges1 = fringe_array(:,fringe_col);
    ridges1 = ridges1';
    ridges1 = ridges1(:)';
    ridges1 = sort(ridges1);
    ridges1 = ridges1(ridges1 ~= 0);
    ridges{n} = ridges1';
    %store only the data on where the center of each fringe is in the cell
    %fringe middles
    fringe_middles{n} = fringe_array(:,6);
    %store the outer bounds (max curvature bounds) in outer bounds
    fringe_col = [4 5];
    bounds1 = fringe_array(:,fringe_col);
    bounds1 = bounds1';
    bounds1 = bounds1(:)';
    bounds1 = sort(bounds1);
    outer_bounds{n} = bounds1';
    %strore the inner bounds (zero curvature/max slope) in inner bounds
    fringe_col = [7 8];
    bounds2 = fringe_array(:,fringe_col);
    bounds2 = bounds2';
    bounds2 = bounds2(:)';
    bounds2 = sort(bounds2);
    inner_bounds{n} = bounds2';
end
%fringe array: [alpha(1) beta(2) a(3) b(4) c(5) d(6) e(7) f(8) g(9) h(10) gamma(11) i(12) j(13) k(14) l(15) m(16)]
%alpha(1) is 2nd derivative peak (left-most peak if not singlet)
%beta(2) is 0 if singlet, second left-most peak if not singlet
%a(3) fringe classification: 0 singlet, 1 doublet, 2 for more than a doublet
%b(4) left max curvature edge
%c(5) right max curvature edge
%d(6) fringe center if not 0
%e(7) left zero curvature edge
%f(8) right zero curvature edge
%g(9) left-most ridge
%h(10) 0 if single ridge, next peak to the right if not
%gamma(11) zero if singlet, right-most neg peak index if not
%i(12) 0 if single ridge, next peak to the right if not
%j(13) 0 if single ridge, next peak to the right if not
%k(14) if more than 2 ridges, the indices of the 1st of the two highest ridges
%l(15) if more than 2 ridges, the indices of the 2nd of the two highest ridges
%m(16) number of ridges within the fringe

%build 2 sets of data: one with all ridges, and one with only fringe middles
a = length(fringes);
fringe_data = [];
ridge_data = [];
middle_data = [];
outer_bounds_data = [];
inner_bounds_data = [];
for n = 1:a
    %make data set with all fringe locations
    fringes1 = fringes2{n};
    c = length(fringes1);
    temp = ones(c,1)*n;
    fringes1 = [temp,fringes1];
    fringe_data = [fringe_data;fringes1];
    %make data set with all ridge locations
    ridges1 = ridges{n};
    c = length(ridges1);
    temp = ones(c,1)*n;
    ridges1 = [temp,ridges1];
    ridge_data = [ridge_data;ridges1];
    %make data set with all fringe middles
    middles1 = fringe_middles{n};
    c = length(middles1);
    temp = ones(c,1)*n;
    middles1 = [temp,middles1];
    middle_data = [middle_data;middles1];
    %make data for all outer bounds
    bounds1 = outer_bounds{n};
    c = length(bounds1);
    temp = ones(c,1)*n;
    bounds1 = [temp,bounds1];
    outer_bounds_data = [outer_bounds_data;bounds1];
    %make data for all inner bounds
    bounds2 = inner_bounds{n};
    c = length(bounds2);
    temp = ones(c,1)*n;
    bounds2 = [temp,bounds2];
    inner_bounds_data = [inner_bounds_data;bounds2];
end
figure;plot(fringe_data(:,1),fringe_data(:,2),'o');title('fringe stuff');
figure;plot(ridge_data(:,1),ridge_data(:,2),'o');title('fringe ridges');
figure;plot(middle_data(:,1),middle_data(:,2),'o');title('fringe middles');
figure;plot(outer_bounds_data(:,1),outer_bounds_data(:,2),'o');title('fringe outer bounds');
figure;plot(inner_bounds_data(:,1),inner_bounds_data(:,2),'o');title('fringe inner bounds');

x=ridge_data(:,1);y=ridge_data(:,2);