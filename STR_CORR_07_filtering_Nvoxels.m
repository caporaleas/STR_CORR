function [AtlasGM_masks_filtered,AtlasWM_masks_filtered,GM_ROIs_flags,WM_ROIs_flags,N_voxels_GM,N_voxels_WM,N_voxel_GM_thr] = STR_CORR_07_filtering_Nvoxels(initpath,Nvoxelfilter,Nvxl_perc,N_GM_ROIs,N_WM_ROIs,V_R1_maps,PVE_masks,AtlasGM_masks,AtlasWM_masks)
% Applying Nvoxel filter (discarding the 20% ROIs with fewer voxels) keeping ROIs with Nvoxels_ROI > 20% 
    if Nvoxelfilter == 1

        cd(initpath)
        folders = dir('pil*');
        N_subj = 20;  % Number of healthy subjects

        disp('Calculating N_voxels for each ROI')
        N_voxels_GM = zeros(N_subj,N_GM_ROIs);
        N_voxels_WM = zeros(N_subj,N_WM_ROIs);

        for ii = 1 : N_subj
            for rr = 1 : N_GM_ROIs
                tmp = V_R1_maps{ii,1}.*PVE_masks{ii,1}.*AtlasGM_masks{rr,1};
                V_R1_masked_arr=tmp(tmp~=0);

                N_voxels_GM(ii,rr) = numel(V_R1_masked_arr);
            end
            for rr = 1 : N_WM_ROIs
                tmp = V_R1_maps{ii,1}.*PVE_masks{ii,2}.*AtlasWM_masks{rr,1};
                V_R1_masked_arr=tmp(tmp~=0);

                N_voxels_WM(ii,rr) = numel(V_R1_masked_arr);
            end
        end

        N_voxels_GM_allsubj = mean(N_voxels_GM,1);
        N_voxel_GM_thr = prctile(N_voxels_GM_allsubj,Nvxl_perc);
        N_voxel_GM_thr = ceil(N_voxel_GM_thr);
        GM_ROIs_flags = N_voxels_GM_allsubj > N_voxel_GM_thr;
        N_ROIs_GM_large_enough = 1 - (N_voxels_GM_allsubj < N_voxel_GM_thr);
        N_GM_ROIs_large_enough = sum(N_ROIs_GM_large_enough);
        sprintf('There are %d GM ROIs large enough out of %d ROIs in the atlas',N_GM_ROIs_large_enough,N_GM_ROIs)

        N_voxels_WM_allsubj = mean(N_voxels_WM,1);
        N_voxel_WM_thr = prctile(N_voxels_WM_allsubj,Nvxl_perc);
        N_voxel_WM_thr = ceil(N_voxel_WM_thr);
        WM_ROIs_flags = N_voxels_WM_allsubj > N_voxel_WM_thr;
        N_ROIs_WM_large_enough = 1 - (N_voxels_WM_allsubj < N_voxel_WM_thr);
        N_WM_ROIs_large_enough = sum(N_ROIs_WM_large_enough);
        sprintf('There are %d WM ROIs large enough out of %d ROIs in the atlas',N_WM_ROIs_large_enough,N_WM_ROIs)

        AtlasGM_masks_filtered = AtlasGM_masks(GM_ROIs_flags);
        AtlasWM_masks_filtered = AtlasWM_masks(WM_ROIs_flags);

        disp('Nvoxels filtering applied')
    else
        AtlasGM_masks_filtered = AtlasGM_masks;
        AtlasWM_masks_filtered = AtlasWM_masks;
        GM_ROIs_flags = ones(180,1);
        disp('No filter on ROI-size applied')
    end
end


