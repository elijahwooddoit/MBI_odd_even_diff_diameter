%algorithm for finding whether a peak is a singlet or a doublet: if
%doublet, there should only be one positive peak between two adjacent
%negative peaks
%classification system: 0=singlet, 1=doublet, 2=skew left, 3=skew right,
%4=hybrid ridge left, 5=hybrid ridge right
%
%
%Classify cases 0, 1, and 6 and fill in doublet center
[~,locs1] = findpeaks(deriv2_ff);
[~,locs2] = findpeaks(deriv2_f);
skip = 0;
eliminate_rows = [];
for n = 1:(length(neg_pks)-1)
    if skip == 1
        skip = 0;
        continue
    elseif skip == 2
        skip = 1;
        continue
    end
    span = neg_pks(n):neg_pks(n+1);
    pos_pks1 = intersect(locs1,span);
    pos_pks2 = intersect(locs2,span);
    if length(pos_pks1) == 1    %code to identify (0&1) singlet vs doublet peaks, also fills in 1,2,3,6 in fringe array
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        fringe_array(n,6) = pos_pks1;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    elseif length(pos_pks2) == 1
        fringe_array(n,2) = neg_pks(n+1);
        fringe_array(n,3) = 1;
        fringe_array(n,6) = pos_pks2;
        eliminate_rows = [eliminate_rows,(n+1)];
        skip = 1;
    end
    if n+2 <= length(neg_pks)   %code to indentify (6) triplet peaks, also fills in 1,2,3,6,9,10,11,12 in fringe array
        span = neg_pks(n+1):neg_pks(n+2);
        pos_pks1 = intersect(locs1,span);
        pos_pks2 = intersect(locs2,span);
        if length(pos_pks1) == 1
            fringe_array(n,3) = 6;
            pks_ind = [neg_pks(n),neg_pks(n+1),neg_pks(n+2)];
            pks_ht = intensity_spectrum(pks_ind);
            low_ind = find(pks_ht == min(pks_ht));
            elim_pk = pks(low_ind);
            fringe_array(n,11) = pks_ind(3);
            fringe_array(n,12) = elim_pk;
            pks = pks_ind;
            pks(low_ind) = [];
            fringe_array(n,6) = round(mean(pks(1),pks(2)));
            fringe_array(n,1) = pks_ind(1);
            fringe_array(n,2) = pks_ind(2);
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
            skip = 2;
            eliminate_rows = [eliminate_rows,(n+1),(n+2)];
            error('!!!triplet detected!!!')
        elseif length(pos_pks2) == 1
            fringe_array(n,3) = 6;
            pks_ind = [neg_pks(n),neg_pks(n+1),neg_pks(n+2)];
            pks_ht = intensity_spectrum(pks_ind);
            low_ind = find(pks_ht == min(pks_ht));
            elim_pk = pks(low_ind);
            fringe_array(n,11) = pks_ind(3);
            fringe_array(n,12) = elim_pk;
            pks = pks_ind;
            pks(low_ind) = [];
            fringe_array(n,6) = round(mean(pks(1),pks(2)));
            fringe_array(n,1) = pks_ind(1);
            fringe_array(n,2) = pks_ind(2);
            fringe_array(n,9) = pks(1);
            fringe_array(n,10) = pks(2);
            skip = 2;
            eliminate_rows = [eliminate_rows,(n+1),(n+2)];
            error('!!!triplet detected!!!')
        end
    end
end
fringe_array(eliminate_rows,:) = [];
%
%classify cases (2&3) skewed peaks, leave second peak empty for these cases, fills in 1,2,3,6 of fringe array
[a,~] = size(fringe_array);
for n = 1:a
    if fringe_array(n,3) == 1
        if deriv2_f(fringe_array(n,6)) < 0
            fringe_array(n,3) = 2;
        end
    end
end
for n = 1:a
    if fringe_array(n,3) == 2
        if intensity_spectrum(fringe_array(n,2)) > intensity_spectrum(fringe_array(n,1))
            fringe_array(n,3) = 3;
            fringe_array(n,9) = 0;
        else
            fringe_array(n,10) = 0;
        end
    end
end
%
%classify cases (4&5) hybrid peaks, define minor peak as the corresponding neg curvature peak, fills in 1,2,3 and minor ridge of fringe array
[~,locs1] = findpeaks(deriv1_ff);
xaxis = zeros(1,length(deriv1_ff));
x = 1:length(xaxis);
[x0,~,~,~] = intersections(x,deriv1_ff,x,xaxis);
locs = round(x0);
for n = 1:a
    if fringe_array(n,3) == 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%need to define 6 in fringe array
        span = fringe_array(n,1):fringe_array(n,2);
        pks1 = intersect(locs1,span);
        zeros0 = intersect(locs,span);
        if length(zeros0) == 1
            fringe_array(n,3) = 4;
            if zeros0 > pks1
                fringe_array(n,3) = 5;
                fringe_array(n,9) = fringe_array(n,1);
            else
                fringe_array(n,10) = fringe_array(n,2);
            end
        end
    end
end




    
    
    
    
    
    
    
    
    
    