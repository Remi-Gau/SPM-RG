% Structural must have same dimension as the functionnal (or the other way
% around if you want smoother images - reslice one or the other);

%% INITIALIZE
close all; clear all; clc;

% Folders for the structural
StructuralFolder = fullfile(pwd, 'Structural');
% Name of the structural file on which the activations will be drawn
StructuralFile = fullfile(StructuralFolder, 'MeanStructural.nii');

ExperimentFolder = fullfile(pwd, 'Experiment');

% Images drawn on a white background
BkGrd = 1;
% Threshold to make the background of the structural white: any voxel with
% a value below the maximum value of the structural divided by this
% threshold will be white.
% You might have to play around a bit with this to get it right.
BkGrdThrshld = 15;

% Enter here the coordinates [X Y Z] of the activations of interested you want
% pictures of. Several activations can entered: one per line.
CoordinatesOfInterest = [-60 10 20];

% if several slices must be plotted around an activation
NumberOfSlices = 1; % How many besides the slice of interest. 
% Must be an odd number : always X slices before and after the slice of interest.

% How many slices should be skipped between the slice of interest amd the 
% previous and following ones. Has to be compatible with the NumberOfSlice
% value
Slice2Skip = 2;  

% Orientation of the slices to be plotted. 1=sagittal, 2=frontal, 3=axial
% Can be any combinations of the 3.
View = 1:3;

% Dimension for the figures: [x y width height] will be adapted later
% depending on volume size
Dimension = repmat([0.3 0.4 0.35 0.35],3,1);

% Number of contrast to plot (between 1 and 3)
NbCon2Plot = 3;
% Actual contrast numbers (these numbers must correspond to the t-map file
% suffix: e.g if you want to plot the values of the file 'spmT_0002.img'
% then put 2 in the vector.
% First contrast will be plotted in red then the next one in green and the
% last in blue
ContrastList = [2 3 14];
% if you want to mask or weight the activations to draw, specify here the 
% mask/weight images corresponding to each contrast.
ROI=cell(3,1);
%   EXAMPLES
%   ROI{1,1} = fullfile(ExperimentFolder, 'PFC_RM_ANOVA_Context_p_0.01.img');                
%   ROI{2,1} = fullfile(ExperimentFolder, 'PFC_RM_ANOVA_Percept_p_0.01.img');
%   ROI{3,1} = fullfile(ExperimentFolder, 'PFC_RM_ANOVA_Interaction_p_0.01.img');

% T Threshold for the peak activation values. Get from SPM results.
Threshold = 2;

% Draw a contour around the activated clusters
Contour = 1;

% Print figures as tif
PrintFig = 0;


               
%% STRUCTURAL
Structural = spm_read_vols(spm_vol(StructuralFile));

% Choose a colormap for the underlay
CM_under = gray(256);
% Set the Min/Max values for structural
AbsMaxStruc = max(abs(Structural(:)));
H_RangeStruc = [0 AbsMaxStruc];
if BkGrd
    Structural(Structural<mean(Structural(:))/BkGrdThrshld) = AbsMaxStruc;
    Color = 'w';
else
    Color = 'k'; %#ok<UNRCH>
end

clear AbsMaxStruc



%% FUNCTIONAL
try
    % load the SPM.mat
    load(fullfile(ExperimentFolder, 'SPM.mat'))
    % display the name of the contrasts to plot
    SPM.xCon(ContrastList(1:NbCon2Plot)).name
catch
end

% Stores the name of the images in cell
% This will use the T images
for ContrastIndex=1:length(ContrastList)
    % Contrast of interest 
    ContrastNumber = ContrastList(ContrastIndex);
    if length(num2str(ContrastNumber))<2
        ContrastNumber = strcat('0', num2str(ContrastNumber));       
    end
    ImgsMat{ContrastIndex,1} = fullfile(ExperimentFolder, strcat('spmT_00', num2str(ContrastNumber), '.img'));
    %ImgsMat{ContrastIndex,1} = fullfile(ExperimentFolder, strcat('con_00', num2str(ContrastNumber), '.img'));
    %ImgsMat{ContrastIndex,1} = fullfile(ExperimentFolder, strcat('Beta_00', num2str(ContrastNumber), '.img'));
end
 
% Get the headers of the images
ImgsInfo = spm_vol(ImgsMat);
ImageDim = ImgsInfo{1,1}.dim;

% Set actual dimension for the figures depending on the dimension of the
% images
Position = Dimension.*[1 1 ImageDim(2)/ImageDim(1) ImageDim(3)/ImageDim(1); ...
                       1 1 1                       ImageDim(3)/ImageDim(1); ...
                       1 1 1                       ImageDim(2)/ImageDim(1)];

% Get the data for each image, weight it
for i=1:size(ImgsInfo,1)
    
    A = spm_read_vols(spm_vol(ImgsMat{i,1}));
    
    if ~isempty(ROI{i,1})
        Mask(:,:,:,i) = spm_read_vols(spm_vol(ROI{i,1}));
        A = Mask(:,:,:,i).*A;
    end

    T_Values(:,:,:,i) = A;
    
    clear A Mask
end

clear ImgsInfo

% Getting the coordinates of the voxel of interest
CoordinatesOfInterest(:,end+1) = 1;
Transf_Matrix = spm_get_space(ImgsMat{1,1});
VoxelsOfInterest = [];
for i=1:size(CoordinatesOfInterest,1) 
    VoxelsOfInterest(i,:) = inv(Transf_Matrix)*CoordinatesOfInterest(i,:)'; %#ok<MINV>
end
clear Transf_Matrix
VoxelsOfInterest(:,4)=[] %#ok<NOPTS>



%% PLOT
% loop does works view per view
FirstFig=1;
for i=1:length(View)
    
    % Figure out which slices to plot
    SliceRange = [];
    if NumberOfSlices==1
        SliceRange = VoxelsOfInterest(:,View(i));
    else
        for j=1:size(VoxelsOfInterest,1)
            SliceRange = [SliceRange ...
                          VoxelsOfInterest(j,View(i))-Slice2Skip*(NumberOfSlices)/2 : ...
                          Slice2Skip : ...
                          VoxelsOfInterest(j,View(i))+Slice2Skip*(NumberOfSlices)/2] %#ok<NOPTS>
        end
    end
    
    % Now slice per slice
    for j=1:length(SliceRange)
        
        SliceNumber = SliceRange(j);
        
        clear T_Map_N_S Final_Overlay
        
        
        % Make a figure: the coordinate mentionned is in voxel space
        switch View(i)
            case 1
                SliceType = 'Sagital_Slice_X=';
                F = figure('name', [SliceType num2str(SliceNumber)], ...
                    'Color', Color, 'Units', 'Normalized', 'Position', Position(View(i),:));

            case 2
                SliceType = 'Frontal_Slice_Y=';
                F = figure('name', [SliceType num2str(SliceNumber)], ...
                    'Color', Color, 'Units', 'Normalized', 'Position', Position(View(i),:));

            case 3
                SliceType = 'Axial_Slice_Z=';
                F = figure('name', [SliceType num2str(SliceNumber)], ...
                    'Color', Color, 'Units', 'Normalized', 'Position', Position(View(i),:));
                               
        end
        axes('Position', [0 0 1 1]); %#ok<LAXES>
        
        
        % Get the values of the structural
        switch View(i)
            case 1
                Underlay = fliplr(rot90(squeeze(Structural(SliceNumber,:,:))));
            case 2
                Underlay = fliplr(rot90(squeeze(Structural(:,SliceNumber,:))));
            case 3
                Underlay = fliplr(rot90(Structural(:,:,SliceNumber)));
        end
        
        % Transform the underlay and beta map to RGB values, based on specified colormaps
        % See function convert_to_RGB() for more information
        U_RGB = convert_to_RGB(Underlay, CM_under, H_RangeStruc);
      
        % Plot the underlay
        layer1 = image(U_RGB);
        hold on;
        
        
        % Get the values of the functional for each contrast
        for k=1:NbCon2Plot
            switch View(i)
                case 1
                    T_Map_N_S(:,:,k) = fliplr(rot90(squeeze(T_Values(SliceNumber,:,:,k))));
                case 2
                    T_Map_N_S(:,:,k) = fliplr(rot90(squeeze(T_Values(:,SliceNumber,:,k))));
                case 3
                    T_Map_N_S(:,:,k) = fliplr(rot90(T_Values(:,:,SliceNumber,k)));
            end
        end
        
        Final_Overlay = repmat(zeros(size(T_Map_N_S(:,:,1))), [1 1 3]);
        for k=1:NbCon2Plot   
            Final_Overlay(:,:,k) = T_Map_N_S(:,:,k)>=Threshold;
        end

        layer2 = image(Final_Overlay); 
        
        % Use the T values to create an alpha map (which must be in [0,1])
        % Right now the alpha map is more an all or none but more elaborate
        % scheme can be developped
        alphamap = logical(sum(Final_Overlay,3));
        set(layer2, 'alphaData', alphamap);
        
        % Add contours
        if Contour
            for k=1:NbCon2Plot
                contour(T_Map_N_S(:,:,k)>=Threshold, 1, 'w', 'linewidth', 2);
            end      
        end
            
        % Other details
        box off
        axis off
        
        % Print the figure
        if PrintFig
            if FirstFig %#ok<UNRCH>
                ExportFormat=hgexport('readstyle','default');
                ExportFormat.Format = 'tiff';
                ExportFormat.Background = Color;
                FirstFig=0;
            end
            hgexport(gcf, strcat(SliceType, num2str(SliceNumber), '.tif'), ExportFormat);
        end
    end
    
end