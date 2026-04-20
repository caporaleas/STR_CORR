function [median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM] = STR_CORR_08_medians(initpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered)

%calculates median of R1 and diff parameter within each ROI and for each
%subj
cd(initpath)
folders = dir('pil*'); 
N_subj = 20;  % Number of healthy subjects

N_par = length(parameters); % Number of metrics

N_GM_ROIs = length(AtlasGM_masks_filtered);
N_WM_ROIs = length(AtlasWM_masks_filtered);

median_R1_GM = zeros(N_subj,N_GM_ROIs);
median_R1_WM = zeros(N_subj,N_WM_ROIs);
median_diff_GM = zeros(N_subj,N_GM_ROIs,N_par-1);
median_diff_WM = zeros(N_subj,N_WM_ROIs,N_par-1);

% GM and WM ROIs
disp('Calculate medians for each ROI')
for ii = 1 : N_subj 
    sprintf('Processing subj %d',ii)

    disp('Calculating R1 medians')
    for rr = 1 : N_GM_ROIs

        tmp = V_R1_maps_filtered{ii,1}.*PVE_masks{ii,1}.*AtlasGM_masks_filtered{rr,1};
        V_R1_masked_arr=tmp(tmp~=0);
        median_R1_GM(ii,rr) = median(V_R1_masked_arr);
        %median_R1_GM(ii,rr) = mean(V_R1_masked_arr);

    end
    for rr = 1 : N_WM_ROIs

        tmp = V_R1_maps_filtered{ii,1}.*PVE_masks{ii,2}.*AtlasWM_masks_filtered{rr,1};
        V_R1_masked_arr=tmp(tmp~=0);
        median_R1_WM(ii,rr) = median(V_R1_masked_arr);

    end

    disp('Calculating diffusion metrics medians')
    for pp = 1 : N_par-1
        for rr = 1 : N_GM_ROIs

            tmp = V_diff_maps_filtered{ii,pp}.*PVE_masks{ii,1}.*AtlasGM_masks_filtered{rr,1};
            V_diff_masked_arr=tmp(tmp~=0);
            median_diff_GM(ii,rr,pp) = median(V_diff_masked_arr);
            %median_diff_GM(ii,rr,pp) = mean(V_diff_masked_arr);
        end

        for rr = 1 : N_WM_ROIs

            tmp = V_diff_maps_filtered{ii,pp}.*PVE_masks{ii,2}.*AtlasWM_masks_filtered{rr,1};
            V_diff_masked_arr=tmp(tmp~=0);
            median_diff_WM(ii,rr,pp) = median(V_diff_masked_arr);

        end
         
    end
end

end