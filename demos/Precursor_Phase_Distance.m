function [distance_dict, x_val_plot] = Precursor_Phase_Distance(all_masks, contiguous_mask, pMap_image, fov, length, pix_wid, csv_name)
% Determine the distance in pixels to count from the edge
dis_dec = (pix_wid / fov) * length;
dis = ceil(dis_dec);
disp('dis');
disp(dis);
% Convert pMap and mask images into usable format
pMap = imread(pMap_image);
% Rescale all intensity values from 0 to 1
R = rescale(pMap, 0, 1);
% Import the general mask
mask_all = imbinarize(imread(all_masks));
% Import the mask outline for edge
I = imread(contiguous_mask);
BW = imbinarize(I);
[row, col] = size(BW);
[BW1, ~] = edge(BW, "sobel");
disp(row)
disp(col)
% Initialize a matrix that is the same size as the original images
results_1 = zeros(row, col);
% Create a loop that only includes values for pixels that are unmasked
for i = 1:row
for j = 1:col
if mask_all(i, j) == 1
results_1(i, j) = R(i, j);
end
end
end
% Create a matrix with all the edge locations for use in the distance counting
edges = [];
for f = 1:row
for g = 1:col
if BW1(f, g) == 1
edges = [edges; f g];
end
end
end
[row_edges, ~] = size(edges);
% Initialize the dictionary to store distance values
distance_dict = containers.Map('KeyType', 'double', 'ValueType', 'any');
for o = 1:dis
distance_dict(o) = [];
end
% Populate the dictionary with individual values
for k = 1:row
if mod(k, 10) == 0
disp(k);
end
for m = 1:col
for p = 1:row_edges
for o = 1:dis
if (k - edges(p, 1)) == o && mask_all(k, m) == 1 && m == edges(p, 2)
distance_dict(o) = [distance_dict(o); results_1(k, m)];
break;
elseif (edges(p, 1) - k) == o && mask_all(k, m) == 1 && m == edges(p, 2)
distance_dict(o) = [distance_dict(o); results_1(k, m)];
break;
elseif (edges(p, 2) - m) == o && mask_all(k, m) == 1 && k == edges(p, 1)
distance_dict(o) = [distance_dict(o); results_1(k, m)];
break;
elseif (m - edges(p, 2)) == o && mask_all(k, m) == 1 && k == edges(p, 1)
distance_dict(o) = [distance_dict(o); results_1(k, m)];
break;
end
end
end
end
end
% Calculate x_val_plot
x_val_plot = (1:dis) * (fov / pix_wid);
disp(x_val_plot)
% Save the dictionary to a CSV file
csv_file_path = fullfile('data/', [csv_name, '.csv']);
fid = fopen(csv_file_path, 'w');
fprintf(fid, 'Distance Bin,Values\n');
for o = 1:dis
values = distance_dict(o);
% Convert the values to a comma-separated string
values_str = strjoin(arrayfun(@num2str, values, 'UniformOutput', false), ',');
fprintf(fid, '%d,%s\n', o, values_str);
end
fclose(fid);
end







