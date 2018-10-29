close all; clear; clc;

ExperimentFolder = fullfile(pwd, 'Non_McGurk_Auditory_Answers');
underlay = fullfile(pwd, 'Structural', 'MeanStructural.nii');
overlay = fullfile(ExperimentFolder, strcat('con_0002.img'));        
alphamap = fullfile(ExperimentFolder, strcat('spmT_0002.img'));  
contour = [];
ROI = [];

contour = fullfile(ExperimentFolder, 'INC_A_sup2_CON_A_p_0.01_c_500.img');
% ROI = fullfile(ExperimentFolder, 'INC_A_sup2_CON_A_p_0.01_c_500.img'); 

opt.Contour = struct();
opt.Contour.Threshold = 1.5;

opt.AlphaMap = struct();
opt.AlphaMap.Threshold = 1;

opt.View = 1:3;

opt.Values2Plot = 0;

opt.TargCoord = [-50 10 25];
opt.NumberOfSlices = 1;
opt.Slice2Skip = 1;


 

spm_up_dual_coding(underlay, overlay, alphamap, contour, ROI, opt)