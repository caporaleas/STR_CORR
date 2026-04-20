function [V_R1_maps_filtered,V_diff_maps_filtered] = STR_CORR_06_filtering_Rsoma(initpath,Rsomafilter,Rsoma_min,V_R1_maps,V_diff_maps,parameters)
%APPLIES Rsoma filtering both in R1 and diffusion maps (eliminates voxels
%where Rsoma is lower than the minimum 
cd(initpath)
folders = dir('pil*');
N_subj = 20;  % Number of healthy subjects

if Rsomafilter == 1
   for ii = 1 : N_subj
       Rsoma = V_diff_maps{ii,5};
       R1tmp = V_R1_maps{ii,1};
       R1tmp(Rsoma < Rsoma_min) = 0;
       V_R1_maps{ii,1} = R1tmp;
       for pp = 1 : length(parameters)-1
           tmp = V_diff_maps{ii,pp};
           tmp(Rsoma < Rsoma_min) = 0;
           V_diff_maps{ii,pp} = tmp;
       end
       
   end
    V_R1_maps_filtered = V_R1_maps;
    V_diff_maps_filtered = V_diff_maps;
    disp('Rsoma filter applied')
else
    V_R1_maps_filtered = V_R1_maps;
    V_diff_maps_filtered = V_diff_maps;
    disp('No filter on Rsoma applied')
end
end
