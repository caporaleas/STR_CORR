function [GM_ROIs,N_GM_ROIs,subGM_ROIs, N_subGM_ROIs,WM_ROIs,N_WM_ROIs,GM_rois_masks,subGM_rois_masks,WM_rois_masks] = STR_CORR_01bis_getting_atlas(img_path_GMatlas,img_path_subGMatlas,img_path_WMatlas,GMROIs_erosion,WMROIs_erosion)
%STR_CORR_01_getting_atlas gets atlas from the paths 

% Loading hcp GM atlas

nii = load_untouch_nii(img_path_GMatlas);
V_hcp_tot = double(nii.img);

GM_ROIs = unique(V_hcp_tot(:));
N_GM_ROIs = numel(GM_ROIs)-1;

% Loading aal2 GM atlas (of which only the subcortical ROIs will be
% retained)
nii = load_untouch_nii(img_path_subGMatlas);
V_aal_tot = double(nii.img);
subGM_labels = [41 42 45 46 75 76 77 78 79 80 81 82]';

for ii = 1 : length(subGM_labels)
    tmp = subGM_labels(ii,1);
    V_aal_tot(V_aal_tot==tmp) = tmp.*1000;
end
V_aal_tot(V_aal_tot<4100) = 0;
V_aal_tot = V_aal_tot./1000;
subGM_ROIs = unique(V_aal_tot(:));
N_subGM_ROIs = numel(subGM_ROIs)-1;

% Loading ICBM WM atlas

nii = load_untouch_nii(img_path_WMatlas);
V_icbm_tot = double(nii.img);

WM_ROIs = unique(V_icbm_tot(:));
N_WM_ROIs = numel(WM_ROIs)-1;

GM_rois_masks = cell(N_GM_ROIs,1);
subGM_rois_masks = cell(N_subGM_ROIs,1);
WM_rois_masks = cell(N_WM_ROIs,1);

for rr = 1 : N_GM_ROIs
    tmp = zeros(size(V_hcp_tot));
    tmp(V_hcp_tot==rr) = 1;
    GM_rois_masks{rr,1} = tmp;
end
for rr = 1 : N_subGM_ROIs
    tmp = zeros(size(V_aal_tot));
    tmp(V_aal_tot==subGM_labels(rr)) = 1;
    subGM_rois_masks{rr,1} = tmp;
end
for rr = 1 : N_WM_ROIs
    tmp = zeros(size(V_icbm_tot));
    tmp(V_icbm_tot==rr) = 1;
    WM_rois_masks{rr,1} = tmp;
end

if GMROIs_erosion == 1
    for rr = 1 : N_GM_ROIs
        tmp = GM_rois_masks{rr,1};
        for ss = 1 : size(tmp,3)
            tmp1 = tmp(:,:,ss);
            bw1 = edge(tmp1);
            tmp1 = tmp1-bw1;
            tmp1(tmp1<0) = 0;
            tmp(:,:,ss) = tmp1;
        end
        GM_rois_masks{rr,1} = tmp;
    end
    disp('GM ROIs erosion applied')

else
    GM_rois_masks{rr,1} = tmp;
    disp('No erosion applied on GM ROIs')

end

if GMROIs_erosion == 1
    for rr = 1 : N_subGM_ROIs
        tmp = subGM_rois_masks{rr,1};
        for ss = 1 : size(tmp,3)
            tmp1 = tmp(:,:,ss);
            bw1 = edge(tmp1);
            tmp1 = tmp1-bw1;
            tmp1(tmp1<0) = 0;
            tmp(:,:,ss) = tmp1;
        end
        subGM_rois_masks{rr,1} = tmp;
    end
    disp('GM ROIs erosion applied on subcortical GM ROIs')

else
    subGM_rois_masks{rr,1} = tmp;
    disp('No erosion applied on subcortical GM ROIs')

end

if WMROIs_erosion == 1
    for rr = 1 : N_WM_ROIs
        tmp = WM_rois_masks{rr,1};
        for ss = 1 : size(tmp,3)
            tmp1 = tmp(:,:,ss);
            bw1 = edge(tmp1);
            tmp1 = tmp1-bw1;
            tmp1(tmp1<0) = 0;
            tmp(:,:,ss) = tmp1;
        end
        WM_rois_masks{rr,1} = tmp;
    end
    disp('WM ROIs erosion applied')

else
    WM_rois_masks{rr,1} = tmp;
    disp('No erosion applied on WM ROIs')

end


end