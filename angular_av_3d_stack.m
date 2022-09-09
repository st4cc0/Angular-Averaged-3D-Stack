clearvars
close all


file_name = dir('*.tif');
tif_file_name = file_name(1).name;
tif_info = imfinfo(tif_file_name);
tif_mass = numel(tif_info);

    %% find center
center_search_idx = round(tif_mass*.2):round(tif_mass*.3);
centers = zeros(length(center_search_idx),2);

n = 1;

for cent_i = center_search_idx
    img_cent = imread(tif_file_name,cent_i);
    img_cent = imadjust(img_cent);
    img_cent_bin = imbinarize(img_cent,0.1);
    img_cent_bin = bwareafilt(img_cent_bin,1,'largest');
    
    region_stats = regionprops(img_cent_bin,'Centroid');
    centers(n,:) = region_stats.Centroid;
    
    n = n + 1;
end

c_x = mean(centers(:,1));
c_y = mean(centers(:,2));
centers = round([c_x,c_y]);
    %% rest
    
for i=1:tif_mass
    img_a = imread(tif_file_name,i);
    for this_phi = -pi:.5:pi
        [img_iso_phi,r_iso_phi] = give_iso_phi(img_a,this_phi,centers);

        [r_iso_phi,r_sort_idx] = sort(r_iso_phi);
        img_iso_phi = img_iso_phi(r_sort_idx);

        r_iso_phi = round(r_iso_phi); 
        [r_iso_unique,id_a,id_b] = unique(r_iso_phi);

        img_iso_unique = zeros(length(r_iso_unique),1);
        for j=1:numel(img_iso_unique)
            av_idx = id_b == j;
            img_iso_unique(j) = mean(img_iso_phi(av_idx));
        end
            %adding lost elements
        if this_phi == -pi
            r_max = max(r_iso_unique);
            r_steady = 0:1:r_max;
            img_z_layer = nan(length(r_steady),1);
            img_z_layer(r_iso_unique) = img_iso_unique;
            1;
        else
            r_max = max(r_iso_unique);

            r_steady = 0:1:r_max;
            img_z_layer_new = nan(length(r_steady),1);
            img_z_layer_new(r_iso_unique) = img_iso_unique;

            img_z_layer = dim_add(img_z_layer_new,img_z_layer);

        end       
    end
    this_img_layer = mean(img_z_layer,2,'omitnan');
    if i == 1
        full_layer = uint8(this_img_layer);
    else
        full_layer = dim_add(uint8(this_img_layer),full_layer);
    end
end

% full_layer(full_layer < 5) = 0;
av_img=flip(full_layer');
imwrite(uint8(av_img),'rot_sym_av.png')
imagesc(av_img)
colormap gray
% shading interp

    


% % figure(2)
% % imshow(img_a)
% % hold on
% % plot(centers(1),centers(2),'xr')

function [img_iso_phi,r_iso_phi] = give_iso_phi(array,this_phi,cent)

    [ni,nj] = size(array);
    err = 0.01;
    img_iso_phi = [];
    r_iso_phi = [];
    
    for i=(1:ni)-cent(2)
        for j=(2:nj)-cent(1)
            if (this_phi-err < atan2(i,j)) && (atan2(i,j) < this_phi+err)
                img_iso_phi = [img_iso_phi, array(i+cent(2),j+cent(1))];
                r_iso_phi = [r_iso_phi, sqrt(i.^2 + j.^2)];
            end
        end
    end
end

function add_up_array = dim_add(new_array,old_array)
    
    [l_a,~]     = size(new_array);
    [l_b,row]   = size(old_array);
        
    if l_a < l_b
        add_up_array = nan(l_b,row+1);
    else
        add_up_array = nan(l_a,row+1);
    end
    
    add_up_array(1:l_b,1:row) = old_array;
    add_up_array(1:l_a,row+1) = new_array;
end