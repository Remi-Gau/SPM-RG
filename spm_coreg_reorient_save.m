function M = spm_coreg_reorient_save(varargin)
% FORMAT M = spm_coreg_reorient_save(TargetScan,SourceScan,OtherImg,flags)
% Nothing fancy just reorient the image one image to another one and saves
% the reorientation matrix in case other images need to be orientated
% later.
%
% TargetScan - handle for reference image (see spm_vol).
% SourceScan - handle for source (moved) image.
% OtherImg - Other images to reorient
% flags - a structure containing the following elements:
%          sep      - optimisation sampling steps (mm)
%                     default: [4 2]
%          params   - starting estimates (6 elements)
%                     default: [0 0 0  0 0 0]
%          cost_fun - cost function string:
%                       'mi'  - Mutual Information
%                       'nmi' - Normalised Mutual Information
%                       'ecc' - Entropy Correlation Coefficient
%                       'ncc' - Normalised Cross Correlation
%                     default: 'nmi'
%          tol      - tolerences for accuracy of each param
%                     default: [0.02 0.02 0.02 0.001 0.001 0.001]
%          fwhm     - smoothing to apply to 256x256 joint histogram
%                     default: [7 7]
%          graphics - display coregistration outputs
%                     default: ~spm('CmdLine')
%
%
% M - reorientation matrix that will move the source image to the target
% image

%%
DateFormat = 'yyyy_mm_dd_HH_MM';

if nargin < 1 || isempty(varargin{1})
    TargetScan = spm_vol(spm_select(1,'image','Select reference image'));
else
    TargetScan = varargin{1};
end

if nargin < 2 || isempty(varargin{2})
    SourceScan = spm_vol(spm_select(Inf,'image','Select moved image(s)'));
else
    SourceScan = varargin{2};
end

if nargin < 3 || isempty(varargin{3})
    ImagesFiles2Process = spm_vol(spm_select(Inf,'image','Select other image(s) to reorient'));
else
    ImagesFiles2Process = varargin{3};
end

if nargin < 4 || isempty(varargin{4})
    flags          = spm_get_defaults('coreg.estimate');
    flags.params   = [0 0 0  0 0 0];
    flags.graphics = ~spm('CmdLine');
else
    flags = varargin{4};
end


%%
fprintf('\n\n')
disp('%%%%%%%%%%%%%%%%')
disp('   COREGISTER   ')
disp('%%%%%%%%%%%%%%%%')
fprintf('\n\n')
% FORMAT x = spm_coreg(TargetScan,SourceScan,flags)
% x     - the parameters describing the rigid body rotation, such that a
%         mapping from voxels in G to voxels in F is attained by:
%         SourceScan.mat\spm_matrix(x(:)')*TargetScan.mat
x = spm_coreg(TargetScan,SourceScan,flags);
M = inv(spm_matrix(x(:)')); %Reorientation matrix


%%
fprintf('\n\n')
disp('%%%%%%%%%%%%%%%%')
disp('    REORIENT    ')
disp('%%%%%%%%%%%%%%%%')
fprintf('\n\n')

spm_progress_bar('Init',size(ImagesFiles2Process,1),...
    'Reading current orientations',...
    'Images loaded');

Mats = zeros(4,4,size(ImagesFiles2Process,1));
for FileInd=1:size(ImagesFiles2Process,1)
    spm_progress_bar('Set',FileInd);
    Mats(:,:,FileInd) = spm_get_space(ImagesFiles2Process{FileInd,1});
end

spm_progress_bar('Clear')

spm_progress_bar('Init',size(ImagesFiles2Process,1),...
    'Writing new orientations',...
    'Images loaded');

spm_get_space(SourceScan, M*spm_get_space(SourceScan)); %#ok<MINV>

for FileInd=1:size(ImagesFiles2Process,1)
    spm_progress_bar('Set',FileInd);
    spm_get_space(ImagesFiles2Process{FileInd,1}, M*Mats(:,:,FileInd)); %#ok<MINV>
end

spm_progress_bar('Clear')

SavedMat = fullfile(pwd, strcat('CoregMatrix_', datestr(now, DateFormat), '.mat'));
fprintf(['Saving coregistration matrix to ' SavedMat '.\n']);
save(SavedMat,'M', 'flags', 'TargetScan', 'SourceScan', 'ImagesFiles2Process');

clear SavedMat ImagesFiles2Process

end
