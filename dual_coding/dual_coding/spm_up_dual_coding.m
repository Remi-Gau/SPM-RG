function spm_up_dual_coding(Underlay, Overlay, Alphamap, Contour, ROI, opt)

% TO DO
% check when only plotting just activation or deactivations
% set default for masking background of underlay

%% check inputs
if isempty(Underlay)
    [spm_dir,~,~] = fileparts(which('spm.m'));
    Underlay = fullfile(spm_dir,'canonical', 'avg305T1.nii'); % TO DO: reslice if not the correct dimension
end
ToCheck{1} = Underlay;

if isempty(Overlay)
    error('I need an overlay image.')
end
ToCheck{2} = Overlay;

if isempty(Alphamap)
    error('I need an alpha map image.')
end
ToCheck{3} = Alphamap;

if ~isempty(Contour)
    opt.Contour.Do=1;
    ToCheck{4} = Contour;
    if isfield(opt.Contour, 'Threshold')
        opt.Contour = rmfield(opt.Contour,'Threshold');
    end
elseif isempty(Contour) && ~isfield(opt.Contour, 'Threshold')
    opt.Contour.Do=0;
end

if ~isempty(ROI)
    ToCheck{5} = ROI;
end

spm_check_orientations(spm_vol(char(ToCheck)))



%% set default options

% Dimension for the figures: [x y width height] will be adapted later
% depending on volume size
Dimension = repmat([0.3 0.4 0.35 0.35],3,1);

% T Threshold for the peak activation values. Get from SPM results.
if ~isfield(opt.AlphaMap, 'Threshold')
    opt.AlphaMap.Threshold = 0;
end
if isfield(opt.Contour, 'Threshold')
    opt.Contour.Do=1;
end

% Plot only the positive (1), negative values (-1) or all of them (0)
if ~isfield(opt, 'Values2Plot')
    opt.Values2Plot = 0;
end
 
% Orientation of the slices to be plotted. 1=sagittal, 2=frontal, 3=axial
% Can be any combinations of the 3.
if ~isfield(opt, 'View')
    opt.View = 1:3;
end

% Enter here the coordinates [X Y Z] of the activations of interested you want
% pictures of. Several activations can entered: one per line.
% if nothing is entered the middle of the volume will be chosen
if ~isfield(opt, 'TargCoord')
    opt.TargCoord = [];
end

% if several slices must be plotted around an activation
if ~isfield(opt, 'NumberOfSlices')
    % How many besides the slice of interest.
    % Must be an odd number : always X slices before and after the slice of interest.
    opt.NumberOfSlices = 3;
    % How many slices should be skipped between the slice of interest amd the
    % previous and following ones. Has to be compatible with the NumberOfSlice
    % value
    opt.Slice2Skip = 2;
end

if ~isfield(opt, 'MaskBkGrd')
    % Images drawn on a white background
    opt.MaskBkGrd.do = 1;

    % Threshold to make the background of the structural white: any voxel with
    % a value below the maximum value of the structural divided by this
    % threshold will be white.
    % You might have to play around a bit with this to get it right.
    opt.MaskBkGrd.Threshold = 15;
end





%% open images
underlay = spm_read_vols(spm_vol(Underlay));
overlay = spm_read_vols(spm_vol(Overlay));
alphamap = spm_read_vols(spm_vol(Alphamap));
if  ~isempty(Contour)
    opt.Contour.Do = 1;
    Contour_Mask = spm_read_vols(spm_vol(Contour));
end
if ~isempty(ROI)
    ROI_Mask = spm_read_vols(spm_vol(ROI));
end

ImgInfo = spm_vol(Overlay);
ImgDim = ImgInfo(1).dim;
opt.Mat = ImgInfo(1).mat;

% Will plot around the center the volume
if isempty(opt.TargCoord)
    opt.TargCoord = round(ImgInfo(1).dim/2);
end
clear ImgInfo




%% Underlay
% Choose a colormap for the underlay
CM_under = gray(256);

% Set the Min/Max values for structural
absmax = max(abs(underlay(:)));
H_RangeUnder = [0 absmax];

% make "background" below a certain value either completely black or white
if opt.MaskBkGrd.do
    underlay(underlay<mean(underlay(:))/opt.MaskBkGrd.Threshold) = absmax;
    Color = 'w'; %#ok<*NASGU>
else
    Color = 'k';
end
clear absmax


%% Overlay                
% Choose a colormap for the overlay
CM_over = seismic(512);

% Remove NaN values and sets them to zero
overlay(isnan(overlay)) = 0; 

if opt.Values2Plot>0
    overlay(overlay<0)=0;
    CM_over(1:floor(size(CM_over,1)/2),:) = [];
elseif opt.Values2Plot<0
    overlay(overlay>0)=0;
    CM_over(ceil(size(CM_over,1)/2):end,:) = [];
end

        
%% Transparency    
% Transparency of the t-value
% Set the Min/Max T-values for alpha coding
A_range = [0 opt.AlphaMap.Threshold];
% Voxels with t-stats of 0 will be completely transparent; 
% Voxels with t-stat magnitudes greater or equal than the threshold will be opaque.
% A_range = [A_range(2)-0.001 A_range(2)];

if opt.Values2Plot>0
    alphamap(alphamap<0)=0;
elseif opt.Values2Plot<0
    alphamap(alphamap>0)=0;
end


%% ROI
if ~isempty(ROI)
    overlay(ROI_Mask==0)=0;
    alphamap(ROI_Mask==0)=0;
    Contour_Mask(ROI_Mask==0)=0;
end


%% Select the slices and decide on the range of values to plot
% Set actual dimension for the figures depending on the dimension of the
% images
Position = Dimension.*[1 1 ImgDim(2)/ImgDim(1) ImgDim(3)/ImgDim(1); ...
                       1 1 1                       ImgDim(3)/ImgDim(1); ...
                       1 1 1                       ImgDim(2)/ImgDim(1)];

[SliceRange, Submax] = slice_range(overlay,opt);
     
% Set the Min/Max values for overlay coding
absmax = max(Submax);
switch opt.Values2Plot
    case 0
        H_Range = [-absmax absmax]; % The colormap is symmetric around zero
    case 1
        H_Range = [0 absmax];
    case -1
        H_Range = [-absmax 0];
end
clear absmax
      

%% Plotting
for i=1:length(opt.View)
    
    SliceRange(opt.View(i),:) 
    
    for j=1:length(SliceRange(opt.View(i),:))
        
        SliceNumber = SliceRange(opt.View(i),j);

        switch opt.View(i)
            case 1
                underlay_N_S = underlay(SliceNumber,:,:);
                overlay_N_S = overlay(SliceNumber,:,:);
                alphamap_N_S = alphamap(SliceNumber,:,:);
                if opt.Contour.Do && ~isempty(Contour)
                    contour_N_S = Contour_Mask(SliceNumber,:,:);
                end
                F = figure('name', ['Sagital_Slice_X=' num2str(SliceNumber)], ...
                    'Color', 'k', 'Units', 'Normalized', 'Position', Position(opt.View(i),:));
                
            case 2
                underlay_N_S = underlay(:,SliceNumber,:);
                overlay_N_S = overlay(:,SliceNumber,:);
                alphamap_N_S = alphamap(:,SliceNumber,:);
                if opt.Contour.Do && ~isempty(Contour)
                    contour_N_S = Contour_Mask(:,SliceNumber,:);
                end
                F = figure('name', ['Sagital_Slice_Y=' num2str(SliceNumber)], ....
                    'Color', 'k', 'Units', 'Normalized', 'Position', Position(opt.View(i),:));
                
            case 3
                underlay_N_S = underlay(:,:,SliceNumber);
                overlay_N_S = overlay(:,:,SliceNumber);
                alphamap_N_S = alphamap(:,:,SliceNumber);
                if opt.Contour.Do && ~isempty(Contour)
                    contour_N_S = Contour_Mask(:,:,SliceNumber);
                end
                F = figure('name', ['Sagital_Slice_Z=' num2str(SliceNumber)], ....
                    'Color', 'k', 'Units', 'Normalized', 'Position', Position(opt.View(i),:));
                
        end
        
        underlay_N_S = fliplr(rot90(squeeze(underlay_N_S)));
        overlay_N_S = fliplr(rot90(squeeze(overlay_N_S)));
        alphamap_N_S = fliplr(rot90(squeeze(alphamap_N_S)));
        if opt.Contour.Do && ~isempty(Contour)
            contour_N_S = fliplr(rot90(squeeze(contour_N_S)));
        elseif opt.Contour.Do
            contour_N_S = alphamap_N_S;
        end
        
        axes('Position', [0 0 1 1]); 
        
        % Transform the underlay and beta map to RGB values, based on specified colormaps
        % See function convert_to_RGB() for more information
        U_RGB = convert_to_RGB(underlay_N_S, CM_under, H_RangeUnder);
        O_RGB = convert_to_RGB(overlay_N_S, CM_over, H_Range);
        
        % Plot the underlay
        layer1 = image(U_RGB); axis image
        hold on;
        % Now, add the Beta difference map as an overlay
        layer2 = image(O_RGB); axis image
        
        % Use the T-statistics to create an alpha map (which must be in [0,1])
        alphamap_N_S = abs(alphamap_N_S);
        alphamap_N_S(alphamap_N_S > A_range(2)) = A_range(2);
        alphamap_N_S(alphamap_N_S < A_range(1)) = 0;
        alphamap_N_S = alphamap_N_S/A_range(2); % Normalize
        % Adjust the alpha values of the overlay
        set(layer2, 'alphaData', alphamap_N_S);
        
        % Add some (black) contours to annotate nominal significance
        hold on
        if opt.Contour.Do 
            if isempty(Contour)
                contour(contour_N_S>=opt.Contour.Threshold, 1, 'k', 'linewidth', 2);
            else
                contour(contour_N_S, 1, 'k', 'linewidth', 2);
            end
        end
        
        % Other details
        box off
        axis off
    end

end


dual_coded_colorbar(A_range, H_Range, CM_over)


end




function IMrgb = convert_to_RGB(IM, cm, cmLIM)
% convert_to_RGB - converts any image to truecolor RGB using a specified colormap  
% USAGE: IMrgb = convert_to_RGB(IM, cm, cmLIM)
% INPUTS: 
%    IM    = the image [m x n]
%    cm    = the colormap [p x 3], optional; default = jet(256)
%    cmLIM = the data limits [min max] to be used in the color-mapping 
%            optional; default = [min(IM) max(IM)]
% OUTPUTS: 
%    IMrgb = the truecolor RGB image [m x n x 3]
% Based on ind2rgb from the Image Processing Toolbox
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------
if nargin < 2, cm = jet(256); end
if nargin < 3, cmLIM = [min(IM(:)) max(IM(:))]; end

IM = IM-cmLIM(1);
IM = IM/(cmLIM(2)-cmLIM(1));
nIND = size(cm,1);
IM = round(IM*(nIND-1));


IM = double(IM)+1;
r = zeros(size(IM)); r(:) = cm(IM,1);
g = zeros(size(IM)); g(:) = cm(IM,2);
b = zeros(size(IM)); b(:) = cm(IM,3);

IMrgb = zeros([size(IM),3]);
% Fill in the r, g, and b channels
IMrgb(:,:,1) = r;
IMrgb(:,:,2) = g;
IMrgb(:,:,3) = b;

end


function dual_coded_colorbar(A_range, H_range, CM_over)
% Create a 2D colorbar for the dual-coded overlay

% Set the labels for the colorbar
hue_label = 'Beta Difference';
alpha_label = '|t|';

%%
G = figure('color', [.5 .5 .5], 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.06, 0.35]);

% x represents the range in alpha (abs(t-stats))
x = linspace(A_range(1), A_range(2), 256);
% y represents the range in hue (beta weight difference)
y = linspace(H_range(1), H_range(2), size(CM_over,1));

[X,Y] = meshgrid(x,y); % Transform into a 2D matrix

axis xy;
box on
xlabel(alpha_label)
ylabel(hue_label)
set(gca, 'Xcolor', 'k', 'Ycolor', 'k')
set(gca, 'YAxisLocation', 'right')

% Plot the colorbar
imagesc(x,y,Y); 
colormap(flipud(CM_over));
alpha(X);
alpha('scaled');

end


function rgb = seismic(n)
% seismic(n) creates a colormap, ranging from dark blue via white to dark red.
%
% Nico Sneeuw
% Munich, 31/08/94

if nargin == 0, n = size(get(gcf,'colormap'),1); end

m = ceil(n/3);
top = ones(m,1);
bot = zeros(m,1);
up = (0:m-1)'/m;
down = flipud(up);

r = [bot; up; 1; top; down];
g = [bot; up; 1; down; bot];
b = [up; top; 1; down; bot];
rgb = [r g b];

% rgb-map has size 4m+1 now. The central part will be extracted.

xlarge = 4*m+1-n;
xblue = round(xlarge/2);
xred = xlarge - xblue;
rgb([1:xblue 4*m-xred+2:4*m+1],:) = [];
end


function [SliceRange, Submax] = slice_range(overlay,opt)

SliceRange = [];
Submax = [];

TargCoord = opt.TargCoord;
Mat = opt.Mat;
View = opt.View;
Slice2Skip = opt.Slice2Skip;
NumberOfSlices = opt.NumberOfSlices;

% Getting the coordinates of the voxel of interest
TargCoord(:,end+1) = 1;
TargVoxel = [];
for i=1:size(opt.TargCoord,1) 
    TargVoxel(i,:) = inv(Mat)*TargCoord(i,:)';  %#ok<*MINV,*AGROW>
end
TargVoxel(:,4)=[];
TargVoxel = round(TargVoxel);

for i=1:length(View)
    
    % Get the range of slices to plot in each dimension
    A = [];
    if NumberOfSlices==1
        SliceRange(View(i),:) = TargVoxel(:,View(i));
    else
        for j=1:size(TargVoxel,1)
            A = [ A ...
                  TargVoxel(j,View(i))-Slice2Skip*(NumberOfSlices) : ...
                  Slice2Skip : ...
                  TargVoxel(j,View(i))+Slice2Skip*(NumberOfSlices)];
        end
        SliceRange(View(i),:) = A;
    end
    
    % Gets the maximum beta value to plot
    switch View(i)
        case 1
           TEMP = abs(overlay(SliceRange(1,:),:,:));
        case 2
           TEMP = abs(overlay(:,SliceRange(2,:),:));
        case 3
           TEMP = abs(overlay(:,:,SliceRange(3,:)));
    end
    Submax = [Submax max(TEMP(:))];
    clear TEMP
    
end

end