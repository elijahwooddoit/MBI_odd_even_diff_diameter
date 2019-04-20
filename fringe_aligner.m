function cent_err_offset = fringe_aligner(Im1,center1,Im2,center2)

%INPUTS
%crop_centers(1)
%cropped_middles{1}
%crop_centers(2)
%cropped middles{2}
%OUTPUTS
%array cent_err_offset: [new center, error, offset from original center]

%center1 = crop_centers(1);
%Im1 = cropped_middles{1};
%center2 = crop_centers(2);
%Im2 = cropped_middles{2};

%calculate the image overlap if the two images are placed with their
%centers on top of one another

%find the overlap error over a span of displacements

%find all the relevant distances and positions and put them in an array
[~,a] = size(Im1);
[~,b] = size(Im2);
distance_matrix = zeros(2,3);
distance_matrix(1,2) = center1;
distance_matrix(2,2) = center2;
distance_matrix(1,1) = center1 - 1;
distance_matrix(2,1) = center2 - 1;
distance_matrix(1,3) = a - center1;
distance_matrix(2,3) = b - center2;
%dist_matrix = [distance from center to left edge, center position, distance from center to right edge]

%create a range of center positions to explore for the second image
dist_matrix2 = distance_matrix;
dist_matrix2(2,:) = 0;
cent_orig = distance_matrix(2,2);
if (cent_orig - 150) > 2
    cent_left = cent_orig - 150;
else
    cent_left = 2;
end
if (cent_orig + 150) < (b-1)
    cent_right = cent_orig + 150;
else
    cent_right = (b-1);
end
cent_range = cent_left:cent_right;

%use loop to measure error at each center position for second image
pos_err = zeros(length(cent_range),2);
%pos_err = [new_center_pos, error]
for n = 1:length(cent_range)
    Image1 = Im1;
    Image2 = Im2;
    crop_matrix = zeros(2,2);
    %crop_matrix: [left crop index, right crop index]
    cent2 = cent_range(n);
    dist_left = cent2 - 1;
    dist_right = b - cent2;
    left_crop = distance_matrix(1,1) - dist_left;
    right_crop = distance_matrix(1,3) - dist_right;
    if left_crop == 0
        %do nothing
    elseif left_crop > 0
        crop_matrix(1,1) = abs(left_crop);
    elseif left_crop < 0
        crop_matrix(2,1) = abs(left_crop);
    end
    if right_crop == 0
        %do nothing
    elseif right_crop > 0
        crop_matrix(1,2) = abs(right_crop);
    elseif right_crop < 0
        crop_matrix(2,2) = abs(right_crop);
    end
    %disp(crop_matrix)
    crop_matrix(:,1) = crop_matrix(:,1) + 1;
    crop_matrix(1,2) = a - crop_matrix(1,2);
    crop_matrix(2,2) = b - crop_matrix(2,2);
    Image1 = Image1(:,crop_matrix(1,1):crop_matrix(1,2));
    Image2 = Image2(:,crop_matrix(2,1):crop_matrix(2,2));
    [c,d] = size(Image1);
    num_of_pix = c*d;
    err = Image1 - Image2;
    err = abs(err);
    err = err.^2;
    err = sum(sum(err));
    err = err / num_of_pix;
    err = err^(1/2);
    pos_err(n,1) = cent2;
    pos_err(n,2) = err;
end
%figure;plot(pos_err(:,1),pos_err(:,2))


[~,idx] = min(pos_err(:,2));
cent_err_offset = pos_err(idx,:);
offst = cent_err_offset(1) - distance_matrix(2,2);
cent_err_offset = [cent_err_offset, offst];

end






