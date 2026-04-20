% Written by Alessandra, Stella Caporale, ITAB 25/10/2024

% This program calculates regional correlation for each subject
% and mean regional correlation (over all subjects) for all regions
% between structural parameters (R1, SANDI metrics, FR from CHARMED model)


% Notes: withouth thresholds on Rsoma, keeping fsoma >= 0.1, percentile >= 20,
% good correlations for fsoma and discrete correlations for Rsoma are
% obtained with CV_GM (Rsoma) >= 0.25

%% clean workspace and setting up the working environment
clc
clearvars
close all

addpath(genpath('/home/stella/Documents/MATLAB/'));  % Add NIfTI toolbox
scriptpath = '/media/nas_rete/Work_stella/Structural_Paper_Scripts/REVISIONS';
addpath(genpath(scriptpath)); % Add other functions/scripts used in the current script
%outpath='/media/nas_rete/Work_ekaterina/SANDI_vs_R1_Glove_mni'; % Set outputpath
outpath='/media/nas_rete/Work_stella/Structural_Paper_CorrMaps/REVISIONS'; % Set outputpath
%outpath='/media/nas_rete/Work_stella/Structural_Paper_CorrMaps/REVISIONS/newR1maps';
initpath = '/media/nas_rete/GLOVE_STUDY/DDC/derivatives';

% list of subjects
done1 = save_subfolders_list_to_txt(initpath,outpath);

% list of parameters to correlate with R1
% parameters={'R1','FR','fneurite','fsoma','fextra','Rsoma','Din','De'};
parameters={'R1','fw','fneurite','fsoma','fextra','Rsoma','Din','De'};
N_par = length(parameters);

%% setting up filters/thresholds
% taking >= threshold
PVEerosion = 1; GMROIs_erosion = 1; WMROIs_erosion = 1;
Rsomafilter = 0; Rsoma_min = 12;
fsomafilter = 1; fsoma_min = 0.1;
Nvoxelfilter = 1; Nvxl_perc = 20;
includeSUBGM = 1; % calculate correlations with the inclusion of subcortical GM ROIs
%CVROIfilter = 1; CVROI_perc = 90;   % 10% of CV is higher than this (in Rsoma map)
%CVsubjfilter = 1; CVsubj_perc = 90;

%% 01) getting atlases

img_path_GMatlas='/storage/shared/Atlas/HCPMMP1_on_MNI152_ICBM2009a_nlin_res_on_MNI_T1_2mm.nii.gz';
img_path_WMatlas='/usr/local/fsl_v6062/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz';
img_path_subGMatlas='/media/nas_rete/GLOVE_STUDY/DDC/scripts/Harvard_Oxford/aal2.nii.gz';

if includeSUBGM == 1
   % do something
   [GM_ROIs,N_GM_ROIs,subGM_ROIs,N_subGM_ROIs,WM_ROIs,N_WM_ROIs,AtlasGM_masks,AtlassubGM_masks,AtlasWM_masks] = STR_CORR_01bis_getting_atlas(img_path_GMatlas,img_path_subGMatlas,img_path_WMatlas,GMROIs_erosion,WMROIs_erosion);
else
    [GM_ROIs,N_GM_ROIs,WM_ROIs,N_WM_ROIs,AtlasGM_masks,AtlasWM_masks] = STR_CORR_01_getting_atlas(img_path_GMatlas,img_path_WMatlas,GMROIs_erosion,WMROIs_erosion);
end

%% 02) getting pve (normalized to MNI) and binarizing them based on the absolute maximum

[PVE_masks] = STR_CORR_02_getting_pve(initpath,outpath,PVEerosion);


%% 03) getting R1 and diffusion metrics

[V_R1_maps, V_diff_maps] = STR_CORR_03_getting_metrics(initpath,parameters);

[V_angle_maps] = STR_CORR_03bis_getting_metrics(initpath);

%[V_R1_maps, V_diff_maps] = STR_CORR_03_getting_metrics_newR1(initpath,parameters);
%/media/nas_rete/Work_davide/Cardiff_MP2RAGE/derivatives/pil002/anat/pil002_desc-R1map_MP2RAGE_brain.nii.gz;

%% 03bis) 
% discarding abnormal R1 values (T1 < 500 ms)
R1_thresh = 2;
for ii = 1 : length(V_R1_maps)
    tmp = V_R1_maps{ii,1};
    tmp(tmp>R1_thresh) = 0;
    V_R1_maps{ii,1} = tmp;
end
disp('Abnormal R1 values have been discarded')
%% 04) showing parametric maps 
% FIGURE

done2 = STR_CORR_04_showing_maps(V_R1_maps,V_diff_maps,outpath);

%% 04bis) showing parametric maps in MNI
% FIGURE
inpath = '/storage/ekaterina/Charmed';
done3 = STR_CORR_04bis_showing_maps(inpath,outpath);

%% 05) applying fsoma filtering (depending on the flags)

[V_R1_maps_filtered,V_diff_maps_filtered] = STR_CORR_05_filtering_fsoma(initpath,fsomafilter,fsoma_min,V_R1_maps,V_diff_maps,parameters);
[V_angle_maps_filtered,V_diff_maps_filtered] = STR_CORR_05_filtering_fsoma(initpath,fsomafilter,fsoma_min,V_angle_maps,V_diff_maps,parameters);

%% 06) applying Rsoma filtering (depending on the flags)

[V_R1_maps_filtered,V_diff_maps_filtered] = STR_CORR_06_filtering_Rsoma(initpath,Rsomafilter,Rsoma_min,V_R1_maps_filtered,V_diff_maps_filtered,parameters);
[V_angle_maps_filtered,V_diff_maps_filtered] = STR_CORR_06_filtering_Rsoma(initpath,Rsomafilter,Rsoma_min,V_angle_maps_filtered,V_diff_maps_filtered,parameters);

%% 07) applying Nvoxels filtering (depending on the flags)
if includeSUBGM == 0
    [AtlasGM_masks_filtered,AtlasWM_masks_filtered,GM_ROIs_flags,WM_ROIs_flags,N_voxels_GM,N_voxels_WM,N_voxels_GM_thr] = STR_CORR_07_filtering_Nvoxels(initpath,Nvoxelfilter,Nvxl_perc,N_GM_ROIs,N_WM_ROIs,V_R1_maps_filtered,PVE_masks,AtlasGM_masks,AtlasWM_masks);
else
    [AtlasGM_masks_filtered,AtlassubGM_masks_filtered,AtlasWM_masks_filtered,subGM_ROIs_flags,WM_ROIs_flags] = STR_CORR_07bis_filtering_Nvoxels(initpath,Nvoxelfilter,Nvxl_perc,N_GM_ROIs,N_subGM_ROIs,N_WM_ROIs,V_R1_maps_filtered,PVE_masks,AtlasGM_masks,AtlassubGM_masks,AtlasWM_masks);
end

WM_ROIs_flags_angle = zeros(1,50);
WM_ROIs_flags_angle(19) = 1;
WM_ROIs_flags_angle(20) = 1;
WM_ROIs_flags_angle(3) = 1;
WM_ROIs_flags_angle(5) = 1;
WM_ROIs_flags_angle(23) = 1;
WM_ROIs_flags_angle(24) = 1;
WM_ROIs_flags_angle(29) = 1;
WM_ROIs_flags_angle(30) = 1;
WM_ROIs_flags_angle = logical(WM_ROIs_flags_angle);
AtlasWM_masks_filtered_angle = AtlasWM_masks(WM_ROIs_flags_angle);

[median_R1_GM,median_R1_WM_angle,median_diff_GM,median_diff_WM_angle] = STR_CORR_08_medians(initpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered_angle);
median_R1_WM_angle_m = median(median_R1_WM_angle,'omitnan');

%% 08) calculating medians for each subject and ROI
if includeSUBGM == 0
    [median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM] = STR_CORR_08_medians(initpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
else
    % 08bis used only pve1 (thalamus is misrepresented)
    %[median_R1_GM,median_R1_subGM,median_R1_WM,median_diff_GM,median_diff_subGM,median_diff_WM,GM_ROI_vol,subGM_ROI_vol,WM_ROI_vol] = STR_CORR_08bis_medians(initpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlassubGM_masks_filtered,AtlasWM_masks_filtered);
    % 08tris uses pve1 and pve2 for subcortical GM ROIs
    [median_R1_GM,median_R1_subGM,median_R1_WM,median_diff_GM,median_diff_subGM,median_diff_WM,GM_ROI_vol,subGM_ROI_vol,WM_ROI_vol] = STR_CORR_08tris_medians(initpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlassubGM_masks_filtered,AtlasWM_masks_filtered);

end

[median_angle_GM,median_angle_subGM,median_angle_WM,median_diff_GM,median_diff_subGM,median_diff_WM,GM_ROI_vol,subGM_ROI_vol,WM_ROI_vol] = STR_CORR_08tris_medians(initpath,parameters,V_angle_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks_filtered,AtlassubGM_masks_filtered,AtlasWM_masks_filtered);

%% 09) calculating correlation coefficients between structural metrics (median +/- median absolute deviation across subjs)
% (?) Do brain regions where R1 has a high median across subjects also have a high median diff metrics?
% FIGURE  
confidence_bounds = 1;
[corr_acrossROIs_GM,corr_acrossROIs_WM] = STR_CORR_09_corr_median_medians_acrossROIs(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath,confidence_bounds);

[corr_acrossROIs_WM_angle] = STR_CORR_09bis_corr_median_medians_acrossROIs(median_R1_WM,median_angle_WM,outpath,confidence_bounds);

%% 10) calculating correlation coefficients between structural metrics (median +/- median absolute deviation across ROIs)
% (?) Do subjects who have a high R1 across ROIs also have a high diff
% metric? (FDR correction)
% FIGURE
confidence_bounds = 1;

if includeSUBGM == 0
    [corr_acrossSubjs_GM,corr_acrossSubjs_WM] = STR_CORR_10_corr_median_medians_acrossSubjs(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath,confidence_bounds);
else
    [corr_acrossSubjs_GM_subGM,corr_acrossSubjs_WM] = STR_CORR_10bis_corr_median_medians_acrossSubjs(median_R1_GM,median_R1_subGM,median_R1_WM,median_diff_GM,median_diff_subGM,median_diff_WM,parameters,outpath,confidence_bounds);
end

%% 11) excluding outliers in diffusion params and recalculating correlation coefficients (FDR correction)
% FIGURE
confidence_bounds = 1;
if includeSUBGM == 0
    [corr_acrossROIs_noout_GM,corr_acrossROIs_noout_WM] = STR_CORR_11_corr_median_medians_acrossROIs_noout(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath,confidence_bounds);
else
    [corr_acrossROIs_noout_GM,corr_acrossROIs_noout_WM] = STR_CORR_11bis_corr_median_medians_acrossROIs_noout(median_R1_GM,median_R1_subGM,median_R1_WM,median_diff_GM,median_diff_subGM,median_diff_WM,parameters,outpath,confidence_bounds);
end

%% EXTRA CONTENT - exploring WM, GM and subGM correlation between fextra and R1
% correlation (excluding outlier per group - GM, subGM, WM
outlierremoval = 0;
[corr_acrossROIs_noout_GM,corr_acrossROIs_noout_WM] = STR_CORR_11tris_corr_median_medians_acrossROIs_noout(outlierremoval,median_R1_GM,median_R1_subGM,median_R1_WM,median_diff_GM,median_diff_subGM,median_diff_WM,parameters,outpath,confidence_bounds);

fextra_subGM = median_diff_subGM(:,:,4);
median_fextra_subGM = median(fextra_subGM,'omitnan');

median_RR1_subGM = median(median_R1_subGM,'omitnan');

fw_WM = median_diff_WM(:,:,1);
median_fw_WM = median(fw_WM,'omitnan');
fextra_WM = median_diff_WM(:,:,4);
median_fextra_WM = median(fextra_WM,'omitnan');
xx = median_fw_WM;
yy = median_fextra_WM;
[r,p] = corrcoef(xx, yy, 'rows','complete');

fw_GM = median_diff_GM(:,:,1);
median_fw_GM = median(fw_GM,'omitnan');
fextra_GM = median_diff_GM(:,:,4);
median_fextra_GM = median(fextra_GM,'omitnan');
xx = median_fw_GM;
yy = median_fextra_GM;
[r,p] = corrcoef(xx, yy, 'rows','complete');

R1_GM = median_R1_GM;
median_R1_GM_m = median(R1_GM,'omitnan');
R1_WM = median_R1_WM;
median_R1_WM_m = median(R1_WM,'omitnan');
xx = median_R1_WM_m;
yy = median_fw_WM;
[r,p] = corrcoef(xx, yy, 'rows','complete');



%% EXTRA CONTENT - correlation bw R1 and GM structural features (ROI volume, curvature, thickness)
% Correlation bw R1 and ROI volume
%{
close all
clc
%CorticalThickness_AcrossSubjects
%GMfeatures_AcrossSubjects
T = readtable('/media/nas_rete/Work_stella/Structural_Paper_Scripts/REVISIONS/HCP_labels.txt', 'Delimiter', ' ', 'ReadVariableNames', false);
GMROI_labelnumbers = table2array(T(:,1));
GMROI_labelnames = table2cell(T(:,2));
GMROI_labelnumbers_filtered = GMROI_labelnumbers(GM_ROIs_flags);
GMROI_labelnames_filtered = GMROI_labelnames(GM_ROIs_flags);
%}
%[corr_acrossROIs_volume_GM,corr_acrossROIs_volume_WM,rho_pcorr,pval_pcorr] = STR_CORR_11_corr_median_medians_acrossROIs_volume(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,GM_ROI_vol,WM_ROI_vol,GM_thickness,GM_MeanR1,GM_GrayVol,GM_MeanCurv,outpath,confidence_bounds);

%% 11bis) excluding outliers in diffusion params based on IQR and recalculating correlation coefficients (FDR correction)
% FIGURE

confidence_bounds = 1;
[corr_acrossROIs_noout2_GM,corr_acrossROIs_noout2_WM] = STR_CORR_11bis_corr_median_medians_acrossROIs_noout_IQR(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath,confidence_bounds);

%% 12) calculating correlation over the entire GM
% FIGURE

confidence_bounds = 1;
[corr_acrossSubjs_allGM,corr_acrossSubjs_allWM,SE_R1_GM] = STR_CORR_12_whole_GM(initpath,outpath,parameters,V_R1_maps_filtered,V_diff_maps_filtered,PVE_masks,AtlasGM_masks,AtlasWM_masks,confidence_bounds);

%% 13) plot historams of r for each diffusion parameter (across subjecs)
% histograms of r for each diffusion parameter (across Subjs)

[corr_for_each_ROI_GM,corr_for_each_ROI_WM] = STR_CORR_13_histograms_acrossSubjs(median_R1_GM,median_R1_WM,median_diff_GM,median_diff_WM,parameters,outpath);
savedvars.corr_ROIs_GM = corr_for_each_ROI_GM;
savedvars.corr_ROIs_WM = corr_for_each_ROI_WM;

close all

cd(outpath);
%% 14) Correlation maps (across subjects)
%{
% A)
[r_map_fn,p_map_fn,rPearson_fn,pPearson_fn] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,2),median_R1_WM,median_diff_WM(:,:,2),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
[r_map_fs,p_map_fs,rPearson_fs,pPearson_fs] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,3),median_R1_WM,median_diff_WM(:,:,3),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
[r_map_fe,p_map_fe,rPearson_fe,pPearson_fe] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,4),median_R1_WM,median_diff_WM(:,:,4),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
[r_map_Rs,p_map_Rs,rPearson_Rs,pPearson_Rs] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,5),median_R1_WM,median_diff_WM(:,:,5),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
[r_map_Di,p_map_Di,rPearson_Di,pPearson_Di] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,6),median_R1_WM,median_diff_WM(:,:,6),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);
[r_map_De,p_map_De,rPearson_De,pPearson_De] = create_r_p_maps_fun(median_R1_GM,median_diff_GM(:,:,7),median_R1_WM,median_diff_WM(:,:,7),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);

% B)
% threshold on p - showing just p<0.05
anat_bkg_path = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz';
h1 = corr_colormaps(anat_bkg_path,r_map_fn,p_map_fn,rPearson_fn,outpath,'fn_corr',230);
k1 = corr_colormaps_pval(anat_bkg_path,r_map_fn,p_map_fn,rPearson_fn,outpath,'fn_corr_pval',230);

% keep on editing this line...
out = imtile({strcat(outpath,'/FigRfn_corr_1.png'),strcat(outpath,'/FigRfn_corr_2.png'),'FigRfn_corr_3.png','FigRfn_corr_4.png','FigRfn_corr_5.png','FigRfn_corr_6.png'},'BackgroundColor','black','BorderSize',[1 1],'GridSize',[1 6]);
h = figure; h.Color = [0 0 0]; 
colormap('jet');
imshow(out,[-1 1]);
c=colorbar;
title(c,'r')
c.Color = [1 1 1]; c.FontSize = 14; c.TickLength = 0;
saveas(h,'FigRfn_allSlices','fig');

out = imtile({'FigPfn_corr_pval_1.png','FigPfn_corr_pval_2.png','FigPfn_corr_pval_3.png','FigPfn_corr_pval_4.png','FigPfn_corr_pval_5.png','FigPfn_corr_pval_6.png'},'BackgroundColor','black','BorderSize',[1 1],'GridSize',[1 6]);
h = figure; h.Color = [0 0 0]; colormap('jet');
imshow(out,[0 0.05]);
c=colorbar;
title(c,'r')
c.Color = [1 1 1]; c.FontSize = 14; c.TickLength = 0;
saveas(h,'FigP_allSlices','fig');

% CHECK across subjects STD
[STD_map_fn,STD_arr_fn] = create_STD_maps_fun(median_R1_GM,median_diff_GM(:,:,2),median_R1_WM,median_diff_WM(:,:,2),AtlasGM_masks,AtlasWM_masks,AtlasGM_masks_filtered,AtlasWM_masks_filtered);

h1 = corr_colormaps_STD(anat_bkg_path,STD_map_fn,STD_arr_fn,outpath,'fn_corr',230);
%}
%% 15) fextra vs R1 - whole brain, whole parenchyma, different tissues
%[PVE_masks_full] = STR_CORR_02_getting_pve(initpath,outpath,0);
% whole brain
h = figure;
for ii = 1 : length(V_R1_maps_filtered)
    subplot(4,5,ii);
    var1 = V_R1_maps_filtered{ii,1};
    var2 = V_diff_maps_filtered{ii,4};
    z = size(var1);
    L = z(1)*z(2)*z(3);
    var1_arr = reshape(var1,L,1);
    var2_arr = reshape(var2,L,1);
    indx = (var1_arr ~= 0) & (var2_arr ~= 0);
    var1_arr_pos = var1_arr(indx);
    var2_arr_pos = var2_arr(indx);
    [r,p] = corrcoef(var1_arr_pos, var2_arr_pos, 'rows','complete');
    plot(var1_arr_pos,var2_arr_pos,'.');
    xlabel('R1 (s^{-1})'); ylabel('fextra');
    hold on
    kk = lsline;
    kk.LineWidth = 2; kk.Color = 'k';
    title(strcat('Subj ',num2str(ii),': r= ',num2str(round(r(2),3,'significant')),'; p= ',num2str(round(p(2),3,'significant'))));   
end
saveas(h,fullfile(outpath,'Whole_brain_fextra_vs_R1'),'fig');

% whole parenchyma
k = figure;
%[PVE_masks_full] = STR_CORR_02_getting_pve(initpath,outpath,0);
% PVE_masks{n,1} = GM, PVE_masks{n,2} = WM, PVE_masks{1,3} = CSF;
for ii = 1 : length(V_R1_maps_filtered)
    
    GMmask = PVE_masks{ii,1};
    WMmask = PVE_masks{ii,2};
    parmask = GMmask+WMmask;
    parmask(parmask==2) = 0;
  
    subplot(4,5,ii);
    var1 = V_R1_maps_filtered{ii,1}.*parmask;
    var2 = V_diff_maps_filtered{ii,4}.*parmask;
    z = size(var1);
    L = z(1)*z(2)*z(3);
    var1_arr = reshape(var1,L,1);
    var2_arr = reshape(var2,L,1);
    indx = (var1_arr ~= 0) & (var2_arr ~= 0);
    var1_arr_pos = var1_arr(indx);
    var2_arr_pos = var2_arr(indx);
    [r,p] = corrcoef(var1_arr_pos, var2_arr_pos, 'rows','complete');
    plot(var1_arr_pos,var2_arr_pos,'g.');
    xlabel('R1 (s^{-1})'); ylabel('fextra');
    hold on
    kk = lsline;
    kk.LineWidth = 2; kk.Color = 'k';
    title(strcat('Subj ',num2str(ii),': r= ',num2str(round(r(2),3,'significant')),'; p= ',num2str(round(p(2),3,'significant'))));   
end
saveas(k,fullfile(outpath,'Whole_parenchyma_fextra_vs_R1'),'fig');

% WM only
k1 = figure;
%[PVE_masks_full] = STR_CORR_02_getting_pve(initpath,outpath,0);
% PVE_masks{n,1} = GM, PVE_masks{n,2} = WM, PVE_masks{1,3} = CSF;
for ii = 1 : length(V_R1_maps_filtered)
    
    WMmask = PVE_masks{ii,2};
   
    subplot(4,5,ii);
    var1 = V_R1_maps_filtered{ii,1}.*WMmask;
    var2 = V_diff_maps_filtered{ii,4}.*WMmask;
    z = size(var1);
    L = z(1)*z(2)*z(3);
    var1_arr = reshape(var1,L,1);
    var2_arr = reshape(var2,L,1);
    indx = (var1_arr ~= 0) & (var2_arr ~= 0);
    var1_arr_pos = var1_arr(indx);
    var2_arr_pos = var2_arr(indx);
    [r,p] = corrcoef(var1_arr_pos, var2_arr_pos, 'rows','complete');
    plot(var1_arr_pos,var2_arr_pos,'r.');
    xlabel('R1 (s^{-1})'); ylabel('fextra');
    hold on
    kk = lsline;
    kk.LineWidth = 2; kk.Color = 'k';
    title(strcat('Subj ',num2str(ii),': r= ',num2str(round(r(2),3,'significant')),'; p= ',num2str(round(p(2),3,'significant'))));   
end
saveas(k1,fullfile(outpath,'Whole_WM_fextra_vs_R1'),'fig');

% GM only
k2 = figure;
%[PVE_masks_full] = STR_CORR_02_getting_pve(initpath,outpath,0);
% PVE_masks{n,1} = GM, PVE_masks{n,2} = WM, PVE_masks{1,3} = CSF;
for ii = 1 : length(V_R1_maps_filtered)
    
    GMmask = PVE_masks{ii,1};
   
    subplot(4,5,ii);
    var1 = V_R1_maps_filtered{ii,1}.*GMmask;
    var2 = V_diff_maps_filtered{ii,4}.*GMmask;
    z = size(var1);
    L = z(1)*z(2)*z(3);
    var1_arr = reshape(var1,L,1);
    var2_arr = reshape(var2,L,1);
    %indx = (var1_arr ~= 0) & (var2_arr ~= 0);
    indx = (var1_arr > 0.6) & (var2_arr ~= 0); % I exclude low R1 values 
    var1_arr_pos = var1_arr(indx);
    var2_arr_pos = var2_arr(indx);
    [r,p] = corrcoef(var1_arr_pos, var2_arr_pos, 'rows','complete');
    plot(var1_arr_pos,var2_arr_pos,'c.');
    xlabel('R1 (s^{-1})'); ylabel('fextra');
    hold on
    kk = lsline;
    kk.LineWidth = 2; kk.Color = 'k';
    title(strcat('Subj ',num2str(ii),': r= ',num2str(round(r(2),3,'significant')),'; p= ',num2str(round(p(2),3,'significant'))));   
end
saveas(k2,fullfile(outpath,'Whole_GM_fextra_vs_R1'),'fig');

%% OLD PRESENTATION OF RESULTS - Correlation coefficients in GM and WM - single FIGURE for WM and GM
%{
corr_for_each_subj_GM = zeros(N_subj,N_par-1);
pvalue_for_each_subj_GM = zeros(N_subj,N_par-1);
corr_for_each_subj_WM = zeros(N_subj,N_par-1);
pvalue_for_each_subj_WM = zeros(N_subj,N_par-1);

for ii = 1 : N_subj
    for pp = 1 : N_par-1
        
        [r,p] = corrcoef(median_R1_GM(ii,:), median_SANDI_GM(ii,:,pp), 'rows','complete');
        
        corr_for_each_subj_GM(ii,pp) = r(2);
        pvalue_for_each_subj_GM(ii,pp) = p(2);

        [r,p] = corrcoef(median_R1_WM(ii,:), median_SANDI_WM(ii,:,pp), 'rows','complete');

        corr_for_each_subj_WM(ii,pp) = r(2);
        pvalue_for_each_subj_WM(ii,pp) = p(2);

    end
end
for ii = 1 : N_subj
    for pp = 1 : N_par-1
        
        [r,p] = corrcoef(median_R1_GM(ii,:), median_SANDI_GM(ii,:,pp), 'rows','complete');
        
        corr_for_each_subj_GM(ii,pp) = r(2);
        pvalue_for_each_subj_GM(ii,pp) = p(2);

        [r,p] = corrcoef(median_R1_WM(ii,:), median_SANDI_WM(ii,:,pp), 'rows','complete');

        corr_for_each_subj_WM(ii,pp) = r(2);
        pvalue_for_each_subj_WM(ii,pp) = p(2);

    end
end

mean_corr_for_each_subj_GM=mean(corr_for_each_subj_GM,1);
std_corr_for_each_subj_GM=std(corr_for_each_subj_GM,1);
SE_corr_for_each_subj_GM=std_corr_for_each_subj_GM./sqrt(N_subj);

mean_corr_for_each_subj_WM=mean(corr_for_each_subj_WM,1);
std_corr_for_each_subj_WM=std(corr_for_each_subj_WM,1);
SE_corr_for_each_subj_WM=std_corr_for_each_subj_WM./sqrt(N_subj);

%{
h = figure;
for pp = 1 : N_par-1
    
    titleplot = strcat(parameter{1,pp+1},' vs R1');
    lst = (1 : 1 : N_subj)';

    subplot(2,7,pp)
    plot(lst,corr_for_each_subj_GM(:,pp),'bo'); title(titleplot);
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_subj_GM(1,pp),'b-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_GM(1,pp)+std_corr_for_each_subj_GM(1,pp)),'b--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_GM(1,pp)-std_corr_for_each_subj_GM(1,pp)),'b--','LineWidth',0.5);
    ylabel('r_{GM}','FontSize',14, 'FontWeight','Bold');
    ylim([-1 1]);
    xlabel('Subjs','FontSize',14,'FontWeight','Bold');
    legend('GM')
    grid on

    subplot(2,7,pp+7)
    plot(lst,corr_for_each_subj_WM(:,pp),'ro');
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_subj_WM(1,pp),'r-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_WM(1,pp)+std_corr_for_each_subj_WM(1,pp)),'r--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_WM(1,pp)-std_corr_for_each_subj_WM(1,pp)),'r--','LineWidth',0.5);
  
    ylabel('r_{WM}','FontSize',14, 'FontWeight','Bold');
    xlabel('Subjs','FontSize',14,'FontWeight','Bold');
    ylim([-1 1]);
    legend('WM')
    grid on
    
    h.Color = [1,1,1];
    
end
saveas(h,fullfile(outpath,'r_Subj_GM_WM'),'fig');
close(h)
%}
% HISTOGRAMS OF r ACROSS SUBJECTS in WM and GM
h = figure;
for pp = 1 : N_par-1
    
    titleplot = strcat(parameters{1,pp+1},' vs R1');
    lst = (1 : 1 : N_subj)';

    subplot(2,7,pp)
    plot(lst,corr_for_each_subj_GM(:,pp),'bo'); title(titleplot);
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_subj_GM(1,pp),'b-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_GM(1,pp)+std_corr_for_each_subj_GM(1,pp)),'b--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_GM(1,pp)-std_corr_for_each_subj_GM(1,pp)),'b--','LineWidth',0.5);
    ylabel('r_{GM}','FontSize',14, 'FontWeight','Bold');
    ylim([-1 1]);
    xlabel('Subjs','FontSize',14,'FontWeight','Bold');
    legend('GM')
    grid on

    subplot(2,7,pp+7)
    plot(lst,corr_for_each_subj_WM(:,pp),'ro');
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_subj_WM(1,pp),'r-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_WM(1,pp)+std_corr_for_each_subj_WM(1,pp)),'r--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_subj_WM(1,pp)-std_corr_for_each_subj_WM(1,pp)),'r--','LineWidth',0.5);
  
    ylabel('r_{WM}','FontSize',14, 'FontWeight','Bold');
    xlabel('Subjs','FontSize',14,'FontWeight','Bold');
    ylim([-1 1]);
    legend('WM')
    grid on
    
    h.Color = [1,1,1];
    
end
saveas(h,fullfile(outpath,'r_Subj_GM_WM'),'fig');
close(h)



%% Region-wise analysis: one point for each GM and WM ROIs - old FIGURE
%median_R1_GM(ii,rr), median_SANDI_GM(ii,rr,pp)
% to edit

corr_for_each_ROI_GM = zeros(N_GM_ROIs,N_par-1);
pvalue_for_each_ROI_GM = zeros(N_GM_ROIs,N_par-1);
corr_for_each_ROI_WM = zeros(N_WM_ROIs,N_par-1);
pvalue_for_each_ROI_WM = zeros(N_WM_ROIs,N_par-1);

for rr = 1 : N_GM_ROIs
    for pp = 1 : N_par-1
        [r,p] = corrcoef(median_R1_GM(:,rr), median_SANDI_GM(:,rr,pp), 'rows','complete');
        
        corr_for_each_ROI_GM(rr,pp) = r(2);
        pvalue_for_each_ROI_GM(rr,pp) = p(2);
    end
end
for rr = 1 : N_WM_ROIs
    for pp = 1 : N_par-1
        [r,p] = corrcoef(median_R1_WM(:,rr), median_SANDI_WM(:,rr,pp), 'rows','complete');

        corr_for_each_ROI_WM(rr,pp) = r(2);
        pvalue_for_each_ROI_WM(rr,pp) = p(2);

    end
end

mean_corr_for_each_ROI_GM=mean(corr_for_each_ROI_GM,1);
std_corr_for_each_ROI_GM=std(corr_for_each_ROI_GM,1);
SE_corr_for_each_ROI_GM=std_corr_for_each_ROI_GM./sqrt(N_GM_ROIs_large_enough);

mean_corr_for_each_ROI_WM=mean(corr_for_each_ROI_WM,1);
std_corr_for_each_ROI_WM=std(corr_for_each_ROI_WM,1);
SE_corr_for_each_ROI_WM=std_corr_for_each_ROI_WM./sqrt(N_WM_ROIs_large_enough);

h = figure;
for pp = 1 : N_par-1
    
    titleplot = strcat(parameters{1,pp+1},' vs R1');
    lst = (1 : 1 : N_GM_ROIs)';

    subplot(2,7,pp)
    plot(lst,corr_for_each_ROI_GM(:,pp),'bo'); title(titleplot);
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_ROI_GM(1,pp),'b-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_ROI_GM(1,pp)+std_corr_for_each_ROI_GM(1,pp)),'b--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_ROI_GM(1,pp)-std_corr_for_each_ROI_GM(1,pp)),'b--','LineWidth',0.5);
    ylabel('r_{GM}','FontSize',14, 'FontWeight','Bold');
    ylim([-1 1]);
    xlabel('GM ROIs','FontSize',14,'FontWeight','Bold');
    legend('GM')
    grid on


    lst = (1 : 1 : N_WM_ROIs)';
    subplot(2,7,pp+7)
    plot(lst,corr_for_each_ROI_WM(:,pp),'ro');
    hold on
    plot(lst,ones(size(lst)).*mean_corr_for_each_ROI_WM(1,pp),'r-','LineWidth',1);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_ROI_WM(1,pp)+std_corr_for_each_ROI_WM(1,pp)),'r--','LineWidth',0.5);
    hold on
    plot(lst,ones(size(lst)).*(mean_corr_for_each_ROI_WM(1,pp)-std_corr_for_each_ROI_WM(1,pp)),'r--','LineWidth',0.5);
  
    ylabel('r_{WM}','FontSize',14, 'FontWeight','Bold');
    xlabel('WM ROIs','FontSize',14,'FontWeight','Bold');
    ylim([-1 1]);
    legend('WM')
    grid on
    
    h.Color = [1,1,1];
    saveas(h,fullfile(outpath,'r_ROIs_GM_WM'),'fig');

end
close(h)

%% SCATTERPLOT y=diffusion par; x=R1, for each ROIs averaged on subjects

N_GM_ROIs_R1 = zeros(N_subj,1);
for ii = 1 : N_subj
    N_GM_ROIs_R1(ii,1) = N_GM_ROIs-sum(isnan(median_R1_GM(ii,:)));
end

N_WM_ROIs_R1 = zeros(N_subj,1);
for ii = 1 : N_subj
    N_WM_ROIs_R1(ii,1) = N_WM_ROIs-sum(isnan(median_R1_WM(ii,:)));
end

N_GM_ROIs_diff = zeros(N_subj,N_par-1);
for pp = 1 : N_par-1
    for ii = 1 : N_subj
        N_GM_ROIs_diff(ii,pp) = N_GM_ROIs-sum(isnan(median_SANDI_GM(ii,:,pp)));
    end
end

N_WM_ROIs_diff = zeros(N_subj,N_par-1);
for pp = 1 : N_par-1
    for ii = 1 : N_subj
        N_WM_ROIs_diff(ii,pp) = N_WM_ROIs-sum(isnan(median_SANDI_WM(ii,:,pp)));
    end
end

plotname = 'Scatterplot_diffpar_vs_R1';

h = figure;
for pp = 1 : N_par-1
    if pp == 5
      median_SANDI_GM_par = median_SANDI_GM(:,:,pp);
      CV_GM_par = std(median_SANDI_GM_par,0,1,'omitnan')./(sqrt(N_GM_ROIs).*mean(median_SANDI_GM_par,1,'omitnan'));
      CV_Rsoma_GM_thr = prctile(CV_GM_par,90);
      [m,n] = find(CV_GM_par > CV_Rsoma_GM_thr);
    end

    subplot(2,4,pp) 

    %mean over all subjects for each region
    mean_GM_R1 = mean(median_R1_GM,1,'omitnan');
    median_SANDI_GM_par = median_SANDI_GM(:,:,pp);
    %mean_diff_GM_par = mean(median_SANDI_GM_par,1,'omitnan');
    
    SE_GM_R1 = std(median_R1_GM,0,1,'omitnan')./sqrt(N_GM_ROIs);
    SE_diff_GM_par = std(median_SANDI_GM_par,0,1,'omitnan')./sqrt(N_GM_ROIs_large_enough);

    median_SANDI_GM_par(:,n) = NaN;
    mean_diff_GM_par = mean(median_SANDI_GM_par,1,'omitnan');

    mean_WM_R1 = mean(median_R1_WM,1,'omitnan');
    median_SANDI_WM_par = median_SANDI_WM(:,:,pp);
    mean_diff_WM_par = mean(median_SANDI_WM_par,1,'omitnan');
    
    SE_WM_R1 = std(median_R1_WM,0,1,'omitnan')./sqrt(N_WM_ROIs);
    SE_diff_WM_par = std(median_SANDI_WM_par,0,1,'omitnan')./sqrt(N_WM_ROIs_large_enough);
   
    [r_GM,p_GM] = corrcoef(mean_GM_R1, mean_diff_GM_par, 'rows','complete');

    corr_coef = round(r_GM(2),3,'significant');
    p_value = round(p_GM(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    [r_WM,p_WM] = corrcoef(mean_WM_R1, mean_diff_WM_par, 'rows','complete');

    corr_coef = round(r_WM(2),3,'significant');
    p_value = round(p_WM(2),3,'significant');

    corr_coef_str_WM = num2str(corr_coef);
    p_value_str_WM = num2str(p_value);

    h.Color = [1 1 1];

    s = errorbar(mean_GM_R1, mean_diff_GM_par, SE_diff_GM_par, SE_diff_GM_par, SE_GM_R1, SE_GM_R1,'s');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0 0.3 0.7];
    
    hold on
    plot(mean_GM_R1,mean_diff_GM_par,'bs','MarkerSize',6);

    %{
    % confidence bounds
    x = mean_R1;
    y = mean_diff_GM_par;
    [p,S] = polyfit(x,y,1);
    xv = linspace(min(x), max(x), 150);
    [y_ext,delta] = polyconf(p,xv,S);
    hold on
    %plot(xv, y_ext, '-b')
    patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none')
    %}

    hold on

    %kk = lsline;
    %kk.Color = 'b'; 
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
    
    hold on

    s = errorbar(mean_WM_R1, mean_diff_WM_par, SE_diff_WM_par, SE_diff_WM_par, SE_WM_R1, SE_WM_R1,'s');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'r';
    s.MarkerFaceColor = [1 0 0];
    
    hold on
    plot(mean_WM_R1,mean_diff_WM_par,'rs','MarkerSize',6);

    hold on

    qq = lsline;
    %qq.Color = {'b','r'}; 
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
    
    %txt = {strcat('r = ',corr_coef_str),strcat('p = ',p_value_str)};
    %annotation('textbox',[.5,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FitBoxToText','on');
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)}, {strcat('r_{WM}=',corr_coef_str_WM,'; p=',p_value_str_WM)});

end
saveas(gcf,fullfile(outpath,plotname),'fig');

close(h)   
%% SCATTERPLOT GM ONLY (TO MAGNIFY THE PLOT) WITH OUTLIERS

plotname = 'Scatterplot_GM_diffpar_vs_R1';

GM_ROIs_labels_num = (1:1:180)';
GM_ROIs_labels_str = num2str(GM_ROIs_labels_num);

h = figure;
for pp = 1 : N_par-1
    
    if pp == 5
      median_SANDI_GM_par = median_SANDI_GM(:,:,pp);
      CV_GM_par = std(median_SANDI_GM_par,0,1,'omitnan')./(sqrt(N_GM_ROIs).*mean(median_SANDI_GM_par,1,'omitnan'));
      CV_Rsoma_GM_thr = prctile(CV_GM_par,90);
      [m,n] = find(CV_GM_par > CV_Rsoma_GM_thr);
    end
    if CVfilter == 1
       median_SANDI_GM_par(:,n) = NaN;
    end

    subplot(2,4,pp) 

    %mean over all subjects for each region
    mean_GM_R1 = mean(median_R1_GM,1,'omitnan');
    
    median_SANDI_GM_par = median_SANDI_GM(:,:,pp);
    outlier_tag = outliers_SANDI(:,pp);
    outlier_tag = outlier_tag';

    mean_diff_GM_par = mean(median_SANDI_GM_par,1,'omitnan');

    SE_GM_R1 = std(median_R1_GM,0,1,'omitnan')./sqrt(N_GM_ROIs);
    SE_diff_GM_par = std(median_SANDI_GM_par,0,1,'omitnan')./sqrt(N_GM_ROIs_large_enough);

    mean_WM_R1 = mean(median_R1_WM,1,'omitnan');
    median_SANDI_WM_par = median_SANDI_WM(:,:,pp);
    mean_diff_WM_par = mean(median_SANDI_WM_par,1,'omitnan');
    
    SE_WM_R1 = std(median_R1_WM,0,1,'omitnan')./sqrt(N_WM_ROIs);
    SE_diff_WM_par = std(median_SANDI_WM_par,0,1,'omitnan')./sqrt(N_WM_ROIs_large_enough);
   
    [r_GM,p_GM] = corrcoef(mean_GM_R1, mean_diff_GM_par, 'rows','complete');

    corr_coef = round(r_GM(2),3,'significant');
    p_value = round(p_GM(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    [r_WM,p_WM] = corrcoef(mean_WM_R1, mean_diff_WM_par, 'rows','complete');

    corr_coef = round(r_WM(2),3,'significant');
    p_value = round(p_WM(2),3,'significant');

    corr_coef_str_WM = num2str(corr_coef);
    p_value_str_WM = num2str(p_value);

    h.Color = [1 1 1];

    s = errorbar(mean_GM_R1, mean_diff_GM_par, SE_diff_GM_par, SE_diff_GM_par, SE_GM_R1, SE_GM_R1,'s');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0 0.3 0.7];
    
    hold on
    plot(mean_GM_R1,mean_diff_GM_par,'bs');


    %{
    % confidence bounds
    x = mean_R1;
    y = mean_diff_GM_par;
    [p,S] = polyfit(x,y,1);
    xv = linspace(min(x), max(x), 150);
    [y_ext,delta] = polyconf(p,xv,S);
    hold on
    %plot(xv, y_ext, '-b')
    patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none')
    %}

    hold on

    kk = lsline;
    kk.Color = 'b'; kk.LineWidth = 2;
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);

    hold on
    mean_diff_GM_par_out = mean_diff_GM_par.*outlier_tag;
    mean_GM_R1_out = mean_GM_R1.*outlier_tag;
    mean_diff_GM_par_out(mean_diff_GM_par_out==0)=NaN;
    plot(mean_GM_R1_out,mean_diff_GM_par_out,'rs');
    text(mean_GM_R1_out, mean_diff_GM_par_out, GM_ROIs_labels_str, 'Vert','bottom', 'Horiz','left', 'FontSize',7)
    
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

end
saveas(gcf,fullfile(outpath,plotname),'fig');

%% SCATTERPLOT GM ONLY (TO MAGNIFY THE PLOT) WITHOUT OUTLIERS

plotname = 'Scatterplot_GM_diffpar_vs_R1_noOutliers';

GM_ROIs_labels_num = (1:1:180)';
GM_ROIs_labels_str = num2str(GM_ROIs_labels_num);

h = figure;
for pp = 1 : N_par-1
    
    subplot(2,4,pp) 

    %mean over all subjects for each region
    mean_GM_R1 = mean(median_R1_GM,1,'omitnan');
    median_SANDI_GM_par = median_SANDI_GM(:,:,pp);
    mean_diff_GM_par = mean(median_SANDI_GM_par,1,'omitnan');

    outlier_tag = outliers_SANDI(:,pp);
    outlier_tag = outlier_tag';
    no_outlier_tag = 1-outlier_tag;

    mean_diff_GM_par = mean_diff_GM_par.*no_outlier_tag;
    mean_GM_R1_out = mean_GM_R1.*no_outlier_tag;
    mean_diff_GM_par(mean_diff_GM_par==0)=NaN;
    N_outlier = sum(isnan(mean_diff_GM_par));
    sprintf('Number of outliers is %f',N_outlier)
    
    SE_GM_R1 = std(median_R1_GM,0,1,'omitnan')./sqrt(N_GM_ROIs);
    SE_diff_GM_par = std(median_SANDI_GM_par,0,1,'omitnan')./sqrt(N_GM_ROIs_large_enough);
    SE_diff_GM_par(mean_diff_GM_par==0)=NaN;
    
    [r_GM,p_GM] = corrcoef(mean_GM_R1, mean_diff_GM_par, 'rows','complete');

    corr_coef = round(r_GM(2),3,'significant');
    p_value = round(p_GM(2),3,'significant');

    corr_coef_str_GM = num2str(corr_coef);
    p_value_str_GM = num2str(p_value);
    
    h.Color = [1 1 1];

    s = errorbar(mean_GM_R1, mean_diff_GM_par, SE_diff_GM_par, SE_diff_GM_par, SE_GM_R1, SE_GM_R1,'s');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0 0.3 0.7];
    
    hold on
    plot(mean_GM_R1,mean_diff_GM_par,'bs');


    %{
    % confidence bounds
    x = mean_R1;
    y = mean_diff_GM_par;
    [p,S] = polyfit(x,y,1);
    xv = linspace(min(x), max(x), 150);
    [y_ext,delta] = polyconf(p,xv,S);
    hold on
    %plot(xv, y_ext, '-b')
    patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], 'b', 'FaceAlpha',0.1, 'EdgeColor','none')
    %}

    hold on

    kk = lsline;
    kk.Color = 'b'; kk.LineWidth = 2;
    ylabel(parameters{1,pp+1},'FontSize',14,'FontWeight','Bold');
    xlabel("R1 (s^{-1})",'FontSize',14,'FontWeight','Bold');
    %xlim([0.5 1]); ylim([0 0.5]);
 
    title({strcat('r_{GM}=',corr_coef_str_GM,'; p=',p_value_str_GM)});

end
saveas(gcf,fullfile(outpath,plotname),'fig');

%}

cd(outpath)