function [r_map_diff,p_map_diff,rPearson_diff,pPearson_diff] = create_r_p_maps_fun(R1_GM,diff_GM,R1_WM,diff_WM,AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered)
% creates regional maps of r and p (Pearson's correlation between Diffusion
% parameters and R1 across subjects)

%only GM
L = length(AtlasGM_masks); 
M = length(AtlasWM_masks);

r_map_diff = zeros(size(AtlasGM_masks{1,1})); 
p_map_diff = zeros(size(AtlasGM_masks{1,1}));

rPearson_diff = zeros(L+M,1); pPearson_diff = zeros(L+M,1);

for rr = 1 : L
    xarr = R1_GM(:,rr);
    yarr = diff_GM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete');
    r_map_diff(AtlasGM_masks_filtered{rr,1}==1) = r;
    p_map_diff(AtlasGM_masks_filtered{rr,1}==1) = p;
    rPearson_diff(rr,1) = r;
    pPearson_diff(rr,1) = p;
end

for rr = 1 : M
    xarr = R1_WM(:,rr);
    yarr = diff_WM(:,rr);
    [r,p] = corr(xarr,yarr,'rows','complete');
    r_map_diff(AtlasWM_masks_filtered{rr,1}==1) = r;
    p_map_diff(AtlasWM_masks_filtered{rr,1}==1) = p;
    rPearson_diff(rr+L,1) = r;
    pPearson_diff(rr+L,1) = p;
    
end
end