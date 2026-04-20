function [STD_map_diff,STD_arr_diff] = create_STD_maps_fun(R1_GM,diff_GM,R1_WM,diff_WM,AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered)
% creates regional maps of r and p (Pearson's correlation between Diffusion
% parameters and R1 across subjects)

%only GM
L = length(AtlasGM_masks); 
M = length(AtlasWM_masks);

STD_map_diff = zeros(size(AtlasGM_masks{1,1})); 
STD_arr_diff = zeros(L+M,1); 

for rr = 1 : L
    yarr = diff_GM(:,rr);
    stddev = std(yarr);
    STD_map_diff(AtlasGM_masks_filtered{rr,1}==1) = stddev;
    STD_arr_diff(rr,1) = stddev;
end

for rr = 1 : M
    yarr = diff_WM(:,rr);
    stddev = std(yarr);
    STD_map_diff(AtlasWM_masks_filtered{rr,1}==1) = stddev; 
    STD_arr_diff(rr+L,1) = stddev;
end
end