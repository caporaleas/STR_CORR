function [median_R1_GM, median_R1_WM, median_diff_GM, median_diff_WM] = STR_CORR_08_medians_diff(initpath, parameters, V_R1_maps_filtered, V_diff_maps_filtered, PVE_masks_all, AtlasGM_masks_filtered, AtlasWM_masks_filtered, N_GM_ROIs_all, N_WM_ROIs_all)

% Calculates median of R1 and diffusion parameters within each ROI and for each subject

cd(initpath)
folders = dir('pil*');
N_subj = 20; % Number of healthy subjects
N_par = length(parameters);  % Number of metrics


N_GM_ROIs_max = 180;
N_WM_ROIs_max = 50;

median_R1_GM = nan(N_subj, N_GM_ROIs_max);
median_R1_WM = nan(N_subj, N_WM_ROIs_max);
median_diff_GM = nan(N_subj, N_GM_ROIs_max, N_par - 1);
median_diff_WM = nan(N_subj, N_WM_ROIs_max, N_par - 1);


disp('Calculate medians for each ROI')
for ii = 1 : N_subj 
    fprintf('Processing subj %d\n', ii)

    AtlasGM_masks_filtered_subj = AtlasGM_masks_filtered{ii,1};
    AtlasWM_masks_filtered_subj = AtlasWM_masks_filtered{ii,1};
    PVE_GMmasks = PVE_masks_all{ii,1};
    PVE_WMmasks = PVE_masks_all{ii,2};

    N_GM_ROIs_subj = N_GM_ROIs_all(ii,1);
    N_WM_ROIs_subj = N_WM_ROIs_all(ii,1);

    if N_GM_ROIs_subj ~= 180
        warning('Subject %d has %d GM ROIs instead of expected %d', ii, N_GM_ROIs_subj, N_GM_ROIs_max);
    end
    if N_WM_ROIs_subj ~= 50
        warning('Subject %d has %d WM ROIs instead of expected %d', ii, N_WM_ROIs_subj, N_WM_ROIs_max);
    end

    disp('Calculating R1 medians')
    for rr = 1 : N_GM_ROIs_subj
        tmp = V_R1_maps_filtered{ii,1} .* PVE_GMmasks .* AtlasGM_masks_filtered_subj{rr,1};
        V_R1_masked_arr = tmp(tmp ~= 0);
        median_R1_GM(ii, rr) = median(V_R1_masked_arr);
    end

    for rr = 1 : N_WM_ROIs_subj
        tmp = V_R1_maps_filtered{ii,1} .* PVE_WMmasks .* AtlasWM_masks_filtered_subj{rr,1};
        V_R1_masked_arr = tmp(tmp ~= 0);
        median_R1_WM(ii, rr) = median(V_R1_masked_arr);
    end

    disp('Calculating diffusion metrics medians')
    for pp = 1 : N_par - 1
        for rr = 1 : N_GM_ROIs_subj
            tmp = V_diff_maps_filtered{ii,pp} .* PVE_GMmasks .* AtlasGM_masks_filtered_subj{rr,1};
            V_diff_masked_arr = tmp(tmp ~= 0);
            median_diff_GM(ii, rr, pp) = median(V_diff_masked_arr);
        end

        for rr = 1 : N_WM_ROIs_subj
            tmp = V_diff_maps_filtered{ii,pp} .* PVE_WMmasks .* AtlasWM_masks_filtered_subj{rr,1};
            V_diff_masked_arr = tmp(tmp ~= 0);
            median_diff_WM(ii, rr, pp) = median(V_diff_masked_arr);
        end
    end
end

end