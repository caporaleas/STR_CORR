function h = corr_colormaps(anat_bkg_path,r_map_GEASL,p_map_GEASL,rPearson_GEASL,outpath,figname,L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

myT1 = load_untouch_nii(anat_bkg_path);
myT1 = double(myT1.img);

% r-maps overlayed on T1 map - AXIAL SLICES saved individually and COLORBAR
r_map_GEASL(isnan(r_map_GEASL)) = 0;
p_map_GEASL(isnan(p_map_GEASL)) = 0;

slices = (20:10:70)';
% GE-ASL
for ss = 1 : length(slices)
    n_slice = slices(ss);
    bkg = rot90(myT1(:,:,n_slice));
    rmap = rot90(r_map_GEASL(:,:,n_slice));
    pmap = rot90(p_map_GEASL(:,:,n_slice)); 

    %all_edges = zeros(size(rmap));
    %{
    for rr = 1 : L
        bw1 = rmap;
        bw1(bw1==rPearson_GEASL(rr,1)) = 1000;
        bw1(bw1==1000) = 1;
        bw1(bw1~=1) = 0;
        se = strel('disk',1);
        bw1ero = imerode(bw1,se);
        roi_edge = bw1-bw1ero;
        all_edges = all_edges + roi_edge;
        all_edges(all_edges>0) = 1;
    end
    %}
    %rmap(pmap>0.05) = 0;
    rmap_in = rmap;
    %rmap_in(all_edges==1) = 0;
    %bkg_edges = all_edges.*max(bkg(:));
    %bkg = bkg.*~bkg_edges;
    
    % this is like zooming in a little bit
    cut=3;   % earlier it was 6
   % bkg = bkg(cut:end-cut,cut:end-cut);
    rmap_in = rmap_in(cut:end-cut,cut:end-cut);
    % end zoom
    h = figure;
    ax1 = axes;
    imagesc(bkg);
    colormap(ax1,'gray');
    ax2 = axes;
    imagesc(ax2,rmap_in,'alphadata',rmap_in~=0);
    colormap(ax2,'jet');
    caxis(ax2,[-1 1]);
    ax1.Visible = 'off';
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    set(gca,'FontSize',20,'FontWeight','bold')
    set(gcf,'color','w');
    %c=colorbar;
    %title(c,'r')
    %c.Color = [1 1 1]; c.FontSize = 14; c.TickLength = 0;
    h.Position = [1 1 800 900];
    h.Color = [0 0 0];
     saveas(h,fullfile(outpath,strcat('FigR',figname,'_',num2str(ss))),'fig');
    exportgraphics(h,fullfile(outpath,strcat('FigR',figname,'_',num2str(ss),'.png')),'Resolution',300)
    if ss == length(slices)
        c = colorbar;
        c.Color = [1 1 1]; c.FontSize = 16; c.TickLength = 0; c.Position = [0.90 0.27 0.0144 0.3];         saveas(h,fullfile(outpath,strcat('FigR',figname,'_colorbar')),'fig');
        exportgraphics(h,fullfile(outpath,strcat('FigR',figname,'_colorbar.png')),'Resolution',300)
         close(h);
    else
         close(h);
    end
end

%% coronal and sagittal slices - R values
% sagittal
n_slice_sag = 45;

bkg_sag = permute(myT1, [2 3 1]);
bkg_sag = rot90(bkg_sag(:,:,n_slice_sag));
rmap_sag = permute(r_map_GEASL, [2 3 1]);
rmap_sag = rot90(rmap_sag(:,:,n_slice_sag));
pmap_sag = permute(p_map_GEASL, [2 3 1]);
pmap_sag = rot90(pmap_sag(:,:,n_slice_sag));
%{
all_edges = zeros(size(rmap_sag));
for rr = 1 : length(rPearson_GEASL)
    bw1 = rmap_sag;
    bw1(bw1==rPearson_GEASL(rr,1)) = 1000;
    bw1(bw1==1000) = 1;
    bw1(bw1~=1) = 0;
    se = strel('disk',1);
    bw1ero = imerode(bw1,se);
    roi_edge = bw1-bw1ero;
    all_edges = all_edges + roi_edge;
    all_edges(all_edges>0) = 1;
end
%}
%rmap_sag(pmap_sag>0.05) = 0;
rmap_in = rmap_sag;

%rmap_in(all_edges==1) = 0;
%bkg_edges = all_edges.*max(bkg_sag(:));
%bkg = bkg_sag.*~bkg_edges;
bkg = bkg_sag;
% this is like zooming in a little bit
cut=3;
%bkg = bkg(cut:end-cut,cut:end-cut);
rmap_in = rmap_in(cut:end-cut,cut:end-cut);
% end zoom
h = figure;
ax1 = axes;
imagesc(bkg);
colormap(ax1,'gray');
ax2 = axes;
imagesc(ax2,rmap_in,'alphadata',rmap_in~=0);
colormap(ax2,'jet');
caxis(ax2,[-1 1]);
ax1.Visible = 'off';
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
%c = colorbar;
%c.Color = [1 1 1]; c.FontSize = 14; c.TickLength = 0;
h.Position = [1 1 1000 800];
h.Color = [0 0 0];
saveas(h,fullfile(outpath,strcat('FigRsag',figname)),'fig');
exportgraphics(h,fullfile(outpath,strcat('FigRsag',figname,'.png')),'Resolution',300)
close(h);


% coronal
n_slice_cor = 55;

bkg_cor = permute(myT1, [3 1 2]);
bkg_cor = rot90(bkg_cor(:,:,n_slice_cor),2);
rmap_cor = permute(r_map_GEASL, [3 1 2]);
rmap_cor = rot90(rmap_cor(:,:,n_slice_cor),2);
pmap_cor = permute(p_map_GEASL, [3 1 2]);
pmap_cor = rot90(pmap_cor(:,:,n_slice_cor),2);
%{
all_edges = zeros(size(rmap_cor));
for rr = 1 : length(rPearson_GEASL)
    bw1 = rmap_cor;
    bw1(bw1==rPearson_GEASL(rr,1)) = 1000;
    bw1(bw1==1000) = 1;
    bw1(bw1~=1) = 0;
    se = strel('disk',1);
    bw1ero = imerode(bw1,se);
    roi_edge = bw1-bw1ero;
    all_edges = all_edges + roi_edge;
    all_edges(all_edges>0) = 1;
end
%}
%rmap_cor(pmap_cor>0.05) = 0;
rmap_in = rmap_cor;
%rmap_in(all_edges==1) = 0;
%rmap_in(pmap_cor > 0.05) = 0;
%bkg_edges = all_edges.*max(bkg_cor(:));
%bkg = bkg_cor.*~bkg_edges;
bkg = bkg_cor;
% this is like zooming in a little bit
cut=3;
%bkg = bkg(cut:end-cut,cut:end-cut);
rmap_in = rmap_in(cut:end-cut,cut:end-cut);
% end zoom
h = figure;
ax1 = axes;
imagesc(bkg);
colormap(ax1,'gray');
ax2 = axes;
imagesc(ax2,rmap_in,'alphadata',rmap_in~=0);
colormap(ax2,'jet');
caxis(ax2,[-1 1]);
ax1.Visible = 'off';
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
%c = colorbar;
%c.Color = [1 1 1]; c.FontSize = 14; c.TickLength = 0;
h.Position = [1 1 1000 800];
h.Color = [0 0 0];
saveas(h,fullfile(outpath,strcat('FigRcor',figname)),'fig');
exportgraphics(h,fullfile(outpath,strcat('FigRcor',figname,'.png')),'Resolution',300)

close(h);

end