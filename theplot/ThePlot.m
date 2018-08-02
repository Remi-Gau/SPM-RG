function ThePlot()

%anat_folder, epi_folder, opt

pdir=pwd;
addpath(pdir);

% if nargin < 2 || isempty(do_td),    do_td = 1;   end;

anat_folder = fullfile('/media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_02/FFX_bu/Structural');
epi_folder = {fullfile('/media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_02/Nifti/NoMoCo/01')};

davg_def = 65;

opt.seg.use = 1;
opt.seg.gm_thres = .5;
opt.seg.wm_thres = .5;
opt.seg.csf_thres = .5;

opt.epi.prefix = 'UR';
opt.epi.prefix = 'UR';


%% Check data availability



% do we have TPMs
if opt.seg.use && isempty(spm_select('FPList',anat_folder,'^c[123].*\.nii$'));
    warning('No tissue probability maps in this folder: %s.\n Did you point to the folder containing the segmented brain.')
    opt.use_seg = 0;
end
% Are the TPMs resliced ?
if opt.seg.use && isempty(spm_select('FPList',anat_folder,'^rc[123].*\.nii$'))
    checkregTPM = 1;
    resliceTPM = 1;
else
    checkregTPM = 0;
    resliceTPM = 0;
end


%% Check that we have tissue probability maps at the EPI resolution
if opt.seg.use
    fprintf('Getting segmentation data')
    
    % Use the mean EPI as reference to check the TPMs otherwise the first
    % image of the first run.
    try
        RefEPI = spm_select('ExtFPList',epi_folder{1},['^mean*' opt.epi.prefix '.*\.nii$'],1);
    catch
        RefEPI = spm_select('ExtFPList',epi_folder{1},['^' opt.epi.prefix '.*\.nii$'],1);
    end
    
    % If they are no resliced TPM we check that the existing ones are
    % coregistered to the EPI
    files = spm_select('FPList',anat_folder,'^c[1].*\.nii$');
    files = cellstr(files);
    files = char(cat(1,RefEPI,files));
    
    if checkregTPM
        fprintf(' Tissue probility maps found but they are not at the EPI resoltuion.')
        fprintf(' Please check that the TPMs are well registered with EPIs.')
        spm_check_registration(spm_vol(files))
    end
    
    % Add user input to decide if to continue or not.
    
    if resliceTPM
        % If OK reslice the TPM
        flags = struct('prefix', 'r' ,'mask', 0, 'mean', 0, 'interp', 0, 'which', 2, 'wrap', [0 0 0]);
        files = spm_select('FPList',anat_folder,'^c\d.*\.nii$');
        files = cellstr(files);
        files = cat(1,RefEPI,files);
        
        spm_reslice(files,flags)
        
        clear files flag
    end
    
    % If they are resliced TPM, we check that they are in the same
    % space/res as the EPIs
    files = spm_select('FPList',anat_folder,'^rc[123].*\.nii$');
    files = cellstr(files);
    files = char(cat(1,RefEPI,files));
    
    if spm_check_orientations(spm_vol(files));
        warning(' The TPM are not in the same space/res as the EPIs.')
        opt.use_seg = 0;
    end
    
    spm_check_registration(spm_vol(files))
    
    clear files
end


%% Get average distance to the surface of the brain
if opt.seg.use
    % Create a surface of the brain using the TPM  if we don't have one
    if isempty(spm_select('FPList',anat_folder,'^rc1.*\.surf\.gii$'))
        files = spm_select('FPList',anat_folder,'^rc[12].*\.nii$');
        spm_surf(files,2);
        clear files
    end
    
    FV = gifti(spm_select('FPList',anat_folder,'^rc1.*\.surf\.gii$'));
    center = FV.vertices(FV.faces(:,:),:);
    center = reshape(center, [size(FV.faces,1) 3 3]);
    center = squeeze(mean(center,2));
    ori_dist = sqrt(sum((center.*-1).^2,2))';
    davg = mean(ori_dist);
    
    clear ori_dist center FV
else
    
    davg = davg_def;
end

fprintf(' Average distance to the cortex surface: %0.2f ', davg)


%% Get data from TPM
if opt.seg.use
    
    gm_seg_hdr = spm_vol(spm_select('FPList',anat_folder,'^rc1.*\.nii$'));
    wm_seg_hdr = spm_vol(spm_select('FPList',anat_folder,'^rc2.*\.nii$'));
    csf_seg_hdr = spm_vol(spm_select('FPList',anat_folder,'^rc3.*\.nii$'));
    
    gm_seg_vol = spm_read_vols(gm_seg_hdr);
    wm_seg_vol = spm_read_vols(wm_seg_hdr);
    csf_seg_vol = spm_read_vols(csf_seg_hdr);
    
    gm_seg_mask = gm_seg_vol>opt.seg.gm_thres;
    prob_gm_mask = gm_seg_vol(gm_seg_mask);
    wm_seg_mask = wm_seg_vol>opt.seg.wm_thres;
    prob_wm_mask = wm_seg_vol(wm_seg_mask);
    csf_seg_mask = csf_seg_vol>opt.seg.csf_thres;
    prob_csf_mask = csf_seg_vol(csf_seg_mask);
    
    sum(sum([gm_seg_mask(:)';wm_seg_mask(:)';csf_seg_mask(:)'])>1)
    
end


%% Get data from motion realignment parameters
cd(epi_folder{1})
img = spm_select('ExtFPList',epi_folder{1},['^' opt.epi.prefix '.*\.nii$'],Inf);
flags = struct('motion_parameters', 'on', 'globals', 'off', 'volume_distance', 'on', 'movie', 'off')
NewRP = spmup_realign_qa(img,flags)

realignment_file = spm_select('FPList',epi_folder{1},'^rp.*\.txt$');

% Frame wise
[FD,RMS] = spmup_FD(realignment_file);
RobustOutlier(FD)


[dxyz]=Pythag_motion(realignment_file);


pr = load(realignment_file);



% total displacement at davg
td = zeros(1,size(pr,1));
for ii = 2:size(pr,1)
    
    dx = pr(ii,1) + (pr(ii,4) .* davg);
    dy = pr(ii,2) + (pr(ii,5) .* davg);
    dz = pr(ii,3) + (pr(ii,6) .* davg);
    td(1,ii) = sqrt(dx.^2 + dy.^2 + dz.^2);
end;


% scan-to-scan discplacement at davg
sts = zeros(1,size(pr,1));
for ii = 2:size(pr,1)
    
    dx = (pr(ii,1)-pr(ii-1,1)) + (pr(ii,4) .* davg -pr(ii-1,4) .* davg);
    dy = (pr(ii,2)-pr(ii-1,2)) + (pr(ii,5) .* davg -pr(ii-1,5) .* davg);
    dz = (pr(ii,3)-pr(ii-1,3)) + (pr(ii,6) .* davg -pr(ii-1,6) .* davg);
    sts(1,ii) = sqrt(dx.^2 + dy.^2 + dz.^2);
end;


% optimize scaling
mot = abs(pr(:,1:3));
deg = abs((pr(:,4:6)*180/pi));
scaleme = ([mot deg td' sts']);
scaleme = [ceil(max(abs(scaleme(:)))) * -1 ceil(max(abs(scaleme(:))))];



















%% Get functional images
cd(epi_folder)

fprintf('Reading functional data')

img = spm_select('ExtFPList',pwd,['^' prefix '.*\.nii$'],Inf);

hdr_epi = spm_vol(img);
vol_epi = spm_read_vols(hdr_epi);













spm_check_orientations(cat(1,hdr_epi,gm_seg_hdr,wm_seg_hdr,csf_seg_hdr));

%%
% flatten and detrend the bold image
dd=size(vol_epi);
im=reshape(vol_epi,[dd(1)*dd(2)*dd(3) dd(4)]);

% mean and trend regressors
r0=ones(1,dd(4));
r1=linspace(0,1,dd(4));

% also set up a temporal mask ignoring the first 4 volumes
tmask=r0;
tmask(1:4)=0;

% remove those terms
im=myregress([r0;r1],im,~~tmask);

% create the basic plot

img = struct('masterflat', im);


gm=img.masterflat(gm_seg_mask(:),:);
[~,I_gm] = sort(prob_gm_mask);
gm=gm(I_gm,:);
wm=img.masterflat(wm_seg_mask(:),:);
[~,I_wm] = sort(prob_wm_mask);
wm=wm(I_wm,:);
csf=img.masterflat(csf_seg_mask(:),:);
[~,I_csf] = sort(prob_csf_mask);
csf=csf(I_csf,:);


sz.gm=size(gm,1);
sz.wm=size(wm,1);
sz.csf=size(csf,1);

redlimz=sz.gm;

allmat=[gm;wm;csf];


sett.gsclr=[0 0 0];
sett.gmclr=[1 .5 0];
sett.wmclr=[1 0 1];
sett.csfclr=[0 0 1];
sett.lw=2;
gmlimz=[-20 20];
xlimz=[1 dd(4)];


%% the plot
close all;
h=figure;
set(h,'position',[100 100 1920 1080]);
set(h,'visible','on');

% plot head motion at top
subplot(5,1,1);
plot(1:length(FD),FD,'r');
xlim([1 length(FD)]);
ylim([0 2]);
set(gca,'ytick',[0 2]);
title('Head motion');
ylabel('Framewise Displacement (mm)');

subplot(5,1,2:5);
imagesc(allmat,gmlimz);
set(gca,'ytick',[]);
colormap(gray);
hh=hline(redlimz(1),'g');
set(hh,'linewidth',2);
xlim(xlimz+.5);
title('Timeseries of in-brain voxels, arranged by SPM segmentation-based compartment masks');
xlabel('Volume #');


ax = gca;
axPos = ax.Position;
axPos(1) = axPos(1)-.01;
axPos(3) = .01;
axes('Position',axPos);

color_bar.gm = repmat(linspace(1,0,sz.gm)', [1,3]).*repmat(sett.gmclr, [sz.gm 1]);
color_bar.wm = repmat(linspace(1,0,sz.wm)', [1,3]).*repmat(sett.wmclr, [sz.wm 1]);;
color_bar.csf = repmat(linspace(1,0,sz.csf)', [1,3]).*repmat(sett.csfclr, [sz.csf 1]);;

color_bar = [color_bar.gm;color_bar.wm;color_bar.csf];
color_bar = reshape(color_bar,[size(color_bar,1),1,3]);
color_bar = repmat(color_bar,[1 50 1]);
image(ones(size(color_bar))-color_bar);

% axis('off')
ylabel('Voxels');
set(gca, 'xtick',[],'ytick',[])


cd(pdir)

return

end

function [resid,pred,b]=myregress(r,tc,varargin)

%  taken from Power's The Plot demo

% presumes vox x time input variable structure

if isempty(varargin)
    % use all timepoints
    b=r'\tc';
    pred=r'*b;
    resid=tc-pred';
    
else
    % use only specified timepoints
    tmask=varargin{1,1};
    b=r(:,tmask)'\tc(:,tmask)';
    pred=r'*b;
    resid=tc-pred';
end
end



function RobustOutlier(x)
% Taken from Cyril Pernet's spmup toolbox
% interquartile range
y=sort(x);
j=floor(length(x)/4 + 5/12);
g=(length(x)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1); % lower quartile
k=length(x)-j+1;
qu=(1-g).*y(k)+g.*y(k-1); % higher quartile
value=qu-ql; % inter-quartile range

% robust outliers
M = median(x);
k=(17.63*length(x)-23.64)/(7.74*length(x)-3.71); % Carling's k
outliers=x<(M-k*value) | x>(M+k*value);

if sum(outliers) > 0
    moutlier_matrix = zeros(size(motion_param,1),sum(outliers));
    indices = find(outliers);
    for i=1:sum(outliers)
        moutlier_matrix(indices(i),i) = 1;
    end
end
end