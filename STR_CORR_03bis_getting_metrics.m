function [V_angle_maps] = STR_CORR_03bis_getting_metrics(initpath)
%STR_CORR_03_getting_metrics gets the V1 from the given path
% Directory with subject data
cd(initpath)
folders = dir('pil*');
N_subj = length(folders);  % Number of subjects
% cell initialization
V_angle_maps = cell(N_subj,1);

% getting V1
for ii = 1 : N_subj  
    
    disp('Filling up cells with V1 metrics')
    disp(ii)

    img_path_V1 = strcat(initpath,'/',folders(ii).name,'/resting/dwi/V12std.nii.gz');
    nii = load_untouch_nii(img_path_V1);
    V_V1_tot = double(nii.img);
    V1x = squeeze((V_V1_tot(:,:,:,1)));
    V1y = squeeze((V_V1_tot(:,:,:,2)));
    V1z = squeeze((V_V1_tot(:,:,:,3)));
    ratio = abs(V1z)./sqrt(V1x.^2 + V1y.^2);
    phirad = atan(ratio);
    phigrad = rad2deg(phirad);
    angle = 90-phigrad;
    V_angle_maps{ii,1} = angle;
        
end

end