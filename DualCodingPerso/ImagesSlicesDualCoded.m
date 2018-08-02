%%
%--------------------------------------------------------------------------
% close all
clear
clc
%--------------------------------------------------------------------------



%% 1. Some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------
% Transparency of the t-value
Transparency = 1; % 1 = yes, 0 = no;
% Use an additional image to plot the contours (e.g an image with only the
% clusters with a FWE-corrected significant size)
Contour = 0;
% Use an additional to mask the Beta and t value to be ploted the contours 
ROI = 1;
% Plot only the positive, negative values or all of them
Values2Plot = 1;

% Threshold for the peak activation values. Get from SPM.
Threshold = 2.42;

% Set the Min/Max T-values for alpha coding
A_range = [0 3];
% Voxels with t-stats of 0 will be completely transparent; 
% voxels with t-stat magnitudes greater or equal than 3.5 will be opaque.
if ~Transparency
    A_range = [A_range(2)-0.01 A_range(2)];
end

% Enter here the coordinates [X Y Z] of the activations of interested you want
% pictures of. Several activations can entered: one per line.
CoordinatesOfInterest = [-50 0 36];

CoordinatesOfInterest(:,end+1) = 1;

% if several slices must be plotted around an activation
NumberOfSlices = 1; % How many besides the slice of interest. Must be an odd 
% number : always X slices before and after the SOI.
Slice2Skip = 1; % How many slices shoudld be skipped between the slide of 
% interest amd the X +/- 1 one.

% Orientation of the slices to be plotted. Can be any combinations of the
% 3.
View = [1 2 3]; % 1 = sagittal, 2 = frontal, 3 = axial

% Contrast of interest.
ContrastNumber = 3;

if length(ContrastNumber)<2
    ContrastNumber = strcat('0', num2str(ContrastNumber));
end

% Variables to find images
Experiment = 'Analysis_Percept_BetweenOnset_100ExpBlocks_TimeDer_200HPF_Despiked';
RFX_Level_Test = '2x2_ANOVA_Auditory_Fused-McGurkinINC_McGurkinCON';

BetaDifferenceFile = fullfile(pwd, 'SecondLevel', Experiment, RFX_Level_Test, ...
                            strcat('con_00', num2str(ContrastNumber), '.img'));                        
                      
TValuesFile = fullfile(pwd, 'SecondLevel', Experiment, RFX_Level_Test, ...
                            strcat('spmT_00', num2str(ContrastNumber), '.img'));                        
                        
StructuralFile = fullfile(pwd, 'SecondLevel', Experiment, RFX_Level_Test, ...
                            'MeanStructuralMasked.nii');
                                          
if Contour==0                  
    ContourFile = TValuesFile;
    else
    ContourFile = fullfile(pwd, 'SecondLevel', Experiment, RFX_Level_Test, ...
                            'test.img');
end
                                    
ROI_File = fullfile(pwd, 'SecondLevel', 'ROI_AAL', 'MNI_Precentral_L.nii');                
                  
                  
SPM_File = fullfile(pwd, 'SecondLevel', Experiment, RFX_Level_Test, 'SPM.mat');           
%--------------------------------------------------------------------------                  



%% Load images and SPM.mat
%--------------------------------------------------------------------------
ImgsMat = char({BetaDifferenceFile ; ...
                TValuesFile ; ...
                StructuralFile; ...
                ContourFile; ...
                ROI_File});
       
ImgsInfo = spm_vol(ImgsMat);

BetaDifference = spm_read_vols(ImgsInfo(1));
BetaDifference(isnan(BetaDifference)) = 0; % Remove NaN values and sets them to zero
if Values2Plot>0
    BetaDifference(BetaDifference<0)=0;
elseif Values2Plot<0
    BetaDifference(BetaDifference>0)=0;
end

TValues = spm_read_vols(ImgsInfo(2));
if Values2Plot>0
    TValues(TValues<0)=0;
elseif Values2Plot<0
    TValues(TValues>0)=0;
end

Structural = spm_read_vols(ImgsInfo(3));

Contour_Mask = spm_read_vols(ImgsInfo(4));

ROI_Mask = spm_read_vols(ImgsInfo(5));
if ROI==1
    BetaDifference(ROI_Mask==0)=0;
    TValues(ROI_Mask==0)=0;
    Contour_Mask(ROI_Mask==0)=0;
end

if (size(BetaDifference)~=size(BetaDifference) ...
        | size(BetaDifference)~=size(TValues) ...
        | size(BetaDifference)~=size(Structural) ...
        | size(BetaDifference)~=size(Contour_Mask) ...
        | size(BetaDifference)~=size(ROI_Mask))
    fprintf('\n\n Images have different sizes. \n\n')
    return
end

load(SPM_File);

% Getting the coordinates of the voxel of interest
Transf_Matrix = spm_get_space(BetaDifferenceFile);
VoxelsOfInterest = [];
for i=1:size(CoordinatesOfInterest,1) 
    VoxelsOfInterest(i,:) = inv(Transf_Matrix)*CoordinatesOfInterest(i,:)';
end
VoxelsOfInterest(:,4)=[]
%--------------------------------------------------------------------------



%% 2. Set some more defaults that will affect the appearance of the image
%--------------------------------------------------------------------------
% Set the Min/Max values for hue coding
absmax = max(abs(BetaDifference(:))); 
H_range = [-absmax absmax]; % The colormap is symmetric around zero

% Set the labels for the colorbar
hue_label = 'Beta Difference (CON_Block - INC_Block)';
alpha_label = '|t|';

% Choose a colormap for the underlay
CM_under = gray(256);

% Choose a colormap for the overlay
CM_over = jet(256);
%--------------------------------------------------------------------------



%% 3. Do the actual plotting
%--------------------------------------------------------------------------
for i=1:length(View)
    
    SliceRange = [];
    if NumberOfSlices==1
        SliceRange = VoxelsOfInterest(:,View(i));
    else
        for j=1:size(VoxelsOfInterest,1)
            SliceRange = [SliceRange ...
                          VoxelsOfInterest(j,View(i))-Slice2Skip*(NumberOfSlices)/2 : ...
                          Slice2Skip : ...
                          VoxelsOfInterest(j,View(i))+Slice2Skip*(NumberOfSlices)/2];
        end
    end
    

    for j=1:length(SliceRange)
        
        SliceNumber = SliceRange(j);
        
        % Gets a single slice of data
        % Bmap_N_S = Difference between betas
        % Tmap_N_S = T-statistics 
        % Pmap_N_S = Binary map 
        % Underlay = Structural image 
         
        switch View(i)
            case 1
                Bmap_N_S = fliplr(rot90(squeeze(BetaDifference(SliceNumber,:,:)))); %'Difference between Novel and Standard betas averaged over 28 subjects'
                Tmap_N_S = fliplr(rot90(squeeze(TValues(SliceNumber,:,:)))); %'T-statistics for the paired t-test comparing Novel and Standard betas'
                Pmap_N_S = fliplr(rot90(squeeze(Contour_Mask(SliceNumber,:,:)))); %'Binary map indicating significance at P<0.001 (fdr corrected)'
                Underlay = fliplr(rot90(squeeze(Structural(SliceNumber,:,:)))); %'Structural image ch2bet from MRIcron, warped to functional data
            case 2
                Bmap_N_S = fliplr(rot90(squeeze(BetaDifference(:,SliceNumber,:))));
                Tmap_N_S = fliplr(rot90(squeeze(TValues(:,SliceNumber,:))));
                Pmap_N_S = fliplr(rot90(squeeze(Contour_Mask(:,SliceNumber,:))));
                Underlay = fliplr(rot90(squeeze(Structural(:,SliceNumber,:))));
            case 3
                Bmap_N_S = fliplr(rot90(BetaDifference(:,:,SliceNumber)));
                Tmap_N_S = fliplr(rot90(TValues(:,:,SliceNumber)));
                Pmap_N_S = fliplr(rot90(Contour_Mask(:,:,SliceNumber)));
                Underlay = fliplr(rot90(Structural(:,:,SliceNumber)));
        end
        
        % Make a figure and set of axes
        F = figure('Color', 'k', 'Units', 'Normalized', 'Position', [0.3, 0.4, 0.2, 0.35]);
        axes('Position', [0 0 1 1]);
        
        % Transform the underlay and beta map to RGB values, based on specified colormaps
        % See function convert_to_RGB() for more information
        U_RGB = convert_to_RGB(Underlay, CM_under);
        O_RGB = convert_to_RGB(Bmap_N_S, CM_over, H_range);
        
        % Plot the underlay
        layer1 = image(U_RGB); axis image
        hold on;
        % Now, add the Beta difference map as an overlay
        layer2 = image(O_RGB); axis image
        
        % Use the T-statistics to create an alpha map (which must be in [0,1])
        alphamap = abs(Tmap_N_S);
        alphamap(alphamap > A_range(2)) = A_range(2);
        alphamap(alphamap < A_range(1)) = 0;
        alphamap = alphamap/A_range(2);
        % Adjust the alpha values of the overlay
        set(layer2, 'alphaData', alphamap);
        
        % Add some (black) contours to annotate nominal significance
        hold on;
        if Contour
            [C, CH] = contour(Pmap_N_S, 1, 'k', 'linewidth', 3);
        else
            [C, CH] = contour(Pmap_N_S, Threshold, 'k', 'linewidth', 3);
        end
        
        % Other details
        box off
        set(gca,'tickdir', 'out');
        % axis off
        
        switch View(i)
            case 1
                SliceType = 'Sagital_Slice_X=';
            case 2
                SliceType = 'Frontal_Slice_Y=';
            case 3
                SliceType = 'Axial_Slice_Z=';
        end
        print(gcf, strcat(SliceType, num2str(SliceNumber), '_', SPM.xCon(1,ContrastNumber).name, '.tif'), '-dtiffnocompression')
        
    end

end
%--------------------------------------------------------------------------



%% 4. Create a 2D colorbar for the dual-coded overlay
%--------------------------------------------------------------------------
G = figure('color', 'w', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.06, 0.35]);
x = linspace(A_range(1), A_range(2), 256);

% x represents the range in alpha (abs(t-stats))
y = linspace(H_range(1), H_range(2), size(CM_over,1));

% y represents the range in hue (beta weight difference)
[X,Y] = meshgrid(x,y); % Transform into a 2D matrix

axis xy;
box on
xlabel(alpha_label)
ylabel(hue_label)
set(gca, 'Xcolor', 'k', 'Ycolor', 'k')
% set(gca, 'YAxisLocation', 'right')

print(gcf, strcat('Activation_', SPM.xCon(1,ContrastNumber).name, '_Legend.eps'), '-dpsc2')

% Plot the colorbar
imagesc(x,y,Y); 
colormap(CM_over);
alpha(X);
alpha('scaled');

print(gcf, strcat('Activation_', SPM.xCon(1,ContrastNumber).name, '_Legend.tif'), '-dtiffnocompression')

axis off

print(gcf, strcat('Activation_', SPM.xCon(1,ContrastNumber).name, '_Legend_NoAxis.tif'), '-dtiffnocompression')

%--------------------------------------------------------------------------


