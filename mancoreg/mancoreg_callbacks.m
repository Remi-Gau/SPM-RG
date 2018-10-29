function mancoreg_callbacks(op,varargin)
%
% Callback routines for mancoreg.m
%
% Change LOG
%
% Version 1.0.1
% Radio button cannot be turned off on matlab for linux (6.5.0). Changed to
% two radio buttons for toggle on/off (12.1.2004, JH)
%
% Version 1.0.2
% Added:    Plot of transformation matrix, values are shown next to sliders
%           and "reset transformation" button (12.1.2004, JH)
%

global st mancoregvar;


% 'move'
% Update the position of the bottom (source) image according to user settings
%----------------------------------------------------------------------------

if strcmp(op,'move'),
        
    angl_pitch=get(mancoregvar.hpitch,'Value');
    angl_roll=get(mancoregvar.hroll,'Value');
    angl_yaw=get(mancoregvar.hyaw,'Value');
    dist_x=get(mancoregvar.hx,'Value');
    dist_y=get(mancoregvar.hy,'Value');
    dist_z=get(mancoregvar.hz,'Value');
    
    set(mancoregvar.hpitch_val,'string',num2str(angl_pitch));
    set(mancoregvar.hroll_val,'string',num2str(angl_roll));
    set(mancoregvar.hyaw_val,'string',num2str(angl_yaw));
    
    set(mancoregvar.hx_val,'string',num2str(dist_x));
    set(mancoregvar.hy_val,'string',num2str(dist_y));
    set(mancoregvar.hz_val,'string',num2str(dist_z));
    
    mancoregvar.sourceimage.premul=spm_matrix([dist_x dist_y dist_z angl_pitch angl_roll angl_yaw 1 1 1 0 0 0]);
    if (get(mancoregvar.htoggle_on,'value')==0) % source is currently displayed
        st.vols{2}.premul=spm_matrix([dist_x dist_y dist_z angl_pitch angl_roll angl_yaw 1 1 1 0 0 0]);
    end
    
    plotmat;
    spm_orthviews('redraw');
    
    return
end


% 'toggle_off'
% Toggles between source and target display in bottom window
%--------------------------------------------------------------------------

if strcmp(op,'toggle_off'),
    
    if (get(mancoregvar.htoggle_off,'value')==0)    % Source is to be displayed
        
        set(mancoregvar.htoggle_off,'value',1);
        
    else
        
        set(mancoregvar.htoggle_on,'value',0);
        st.vols{2}=mancoregvar.sourceimage;
        spm_orthviews('redraw');
        
        
    end
    
    return
end


% 'toggle_on'
% Toggles between source and target display in bottom window
%--------------------------------------------------------------------------

if strcmp(op,'toggle_on'),
    
    if (get(mancoregvar.htoggle_on,'value')==0)    % Source is to be displayed
        
        set(mancoregvar.htoggle_on,'value',1);
        
    else
        set(mancoregvar.htoggle_off,'value',0);
        mancoregvar.sourceimage=st.vols{2}; % Backup current state
        st.vols{2}=st.vols{1};
        st.vols{2}.ax=mancoregvar.sourceimage.ax;   % These have to stay the same
        st.vols{2}.window=mancoregvar.sourceimage.window;
        st.vols{2}.area=mancoregvar.sourceimage.area;
        spm_orthviews('redraw');
    end
    
    return
end

% 'reset'
% Resets transformation matrix
%--------------------------------------------------------------------------

if strcmp(op,'reset'),
    
    set(mancoregvar.hpitch,'Value',0);
    set(mancoregvar.hroll,'Value',0);
    set(mancoregvar.hyaw,'Value',0);
    
    set(mancoregvar.hx,'Value',0);
    set(mancoregvar.hy,'Value',0);
    set(mancoregvar.hz,'Value',0);
    
    set(mancoregvar.hpitch_val,'string','0');
    set(mancoregvar.hroll_val,'string','0');
    set(mancoregvar.hyaw_val,'string','0');
    
    set(mancoregvar.hx_val,'string','0');
    set(mancoregvar.hy_val,'string','0');
    set(mancoregvar.hz_val,'string','0');
    
    mancoregvar.sourceimage.premul=spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]);
    if (get(mancoregvar.htoggle_on,'value')==0) % source is currently displayed
        st.vols{2}.premul=spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]);
    end
    
    plotmat;
    spm_orthviews('redraw');
    
    
    return
end


% 'apply'
% Apply transformation to a selected set of images
%--------------------------------------------------------------------------

if strcmp(op,'apply'),
    
    mat = getmat(mancoregvar);
    
    spm_defaults;
    mat = spm_matrix(mat);
    
    % The following is copied from spm_image.m
    if det(mat)<=0
        spm('alert!','This will flip the images',mfilename,0,1);
    end;
    P = spm_select(Inf, 'image','Images to reorient');
    
    
    DateFormat = 'yyyy_mm_dd_HH_MM';
    M=mat; %#ok<NASGU>
    SavedMat = fullfile(pwd, strcat('ReorientMatrix_', datestr(now, DateFormat), '.mat'));
    fprintf(['Saving reorient matrice to ' SavedMat '.\n']);
    save(SavedMat,'M', 'P');
    clear M DateFormat SavedMat
    
    
    
    %fprintf('Skipping image selection! Preselected images ...\n');
    %load d:\mrdata\functional\prediction\v1_chris\P.mat
    
    Mats = zeros(4,4,size(P,1));
    
    for i=1:size(P,1)
        tmp = sprintf('Reading current orientations... %.0f%%.\n', i/size(P,1)*100 );
        fprintf('%s',tmp)
        
        Mats(:,:,i) = spm_get_space(P(i,:));
        spm_progress_bar('Set',i);
        
        fprintf('%s',char(sign(tmp)*8))
    end;
    
    
    for i=1:size(P,1)
        tmp = sprintf('Reorienting images... %.0f%%.\n', i/size(P,1)*100 );
        fprintf('%s',tmp)
        
        spm_get_space(P(i,:),mat*Mats(:,:,i));
        
        fprintf('%s',char(sign(tmp)*8))
    end;
    
    
    tmp = spm_get_space([st.vols{1}.fname ',' num2str(st.vols{1}.n)]);
    if sum((tmp(:)-st.vols{1}.mat(:)).^2) > 1e-8,
        spm_image('init',st.vols{1}.fname);
    end;
    return;
end;


% 'plotmat'
% Display joint histograms
%--------------------------------------------------------------------------

if strcmp(op,'dispcoreghist'),
    
    [VG,VF] = getimg(mancoregvar);
    
    x = getmat(mancoregvar);
    
    coreghist = spm_figure('Findwin','coreghist');
    if isempty(coreghist)
        coreghist=spm_figure('Create','coreghist','Coregistration histogram','on');
        if isempty(coreghist),
            error('Cant create graphics window');
        end;
    else
        spm_figure('Clear','coreghist');
    end;

    ax  = axes('Position',[0.1 0.7 0.35 0.3],'Visible','off','Parent',coreghist);

    H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\VG.mat,[1 1 1]);
    tmp = log(H+1);
    
    image(tmp*(64/max(tmp(:))),'Parent',ax');
    
    set(ax,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
        'XTick',[],'YTick',[]);
    title('Original Joint Histogram','Parent',ax);
    xlabel(spm_file(mancoregvar.targetimage.fname,'short22'),'Parent',ax,'Interpreter','none');
    ylabel(spm_file(mancoregvar.sourceimage.fname,'short22'),'Parent',ax,'Interpreter','none');
    
    
    ax  = axes('Position',[0.6 0.7 0.35 0.3],'Visible','off','Parent',coreghist);
    
    H   = spm_hist2(VG.uint8,VF.uint8,VF.mat\spm_matrix(x(:)')*VG.mat,[1 1 1]);   
    tmp = log(H+1);
    
    image(tmp*(64/max(tmp(:))),'Parent',ax');
    
    set(ax,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal',...
        'XTick',[],'YTick',[]);
    title('Final Joint Histogram','Parent',ax);
    xlabel(spm_file(mancoregvar.targetimage.fname,'short22'),'Parent',ax,'Interpreter','none');
    ylabel(spm_file(mancoregvar.sourceimage.fname,'short22'),'Parent',ax,'Interpreter','none');
    
    return
end


% 'costfunc'
% computes the different type of cost function
%--------------------------------------------------------------------------

if strcmp(op,'costfunc'),
    
    cfname ={
        'Mutual Information'; ...
        'Entropy Correlation Coefficient';...
        'Normalised Mutual Information';...
        'Normalised Cross Correlation'};
    cf = {'mi';'ecc';'nmi';'ncc'};

    [VG,VF] = getimg(mancoregvar);
    
    x = getmat(mancoregvar);
    
    if ~isfield('mancoregvar', 'cf')
        for i=1:numel(cf)
            mancoregvar.cf(i,1) = optfun(1,VG,VF,[1 1 1],cf{i},[7 7]);
        end
         mancoregvar.mat{1,1} = spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]);
    end
    
    Iter = size(mancoregvar.cf,2)+1;
    for i=1:numel(cf)
        mancoregvar.cf(i,Iter) = optfun(x,VG,VF,[1 1 1],cf{i},[7 7]);
    end
    mancoregvar.mat{1,Iter} = x;
    
    
%     cfhist = spm_figure('Findwin','cfhist');
%     if isempty(cfhist)
        cfhist=spm_figure('Create','cfhist','Cost functions','on');
%         if isempty(cfhist),
%             error('Cant create graphics window');
%         end;
%     else
%         spm_figure('Clear','cfhist');
%     end;

    for i=1:numel(cf)
        ax  = axes('Position',[0.1 .15+(i-1)*0.2 0.8 0.15],'Visible','off','Parent',cfhist);
        plot(mancoregvar.cf(i,:))
        set(ax,'XTick',1:Iter,'XTickLabel',1:Iter);
        title(cfname{i},'Parent',ax);
    end

    return
end


% 'plotmat'
% Plot matrix notation of transformation
%--------------------------------------------------------------------------

if strcmp(op,'plotmat'),
    
    plotmat;
    return
end

% None of the op strings matches

fprintf('WARNING: mancoreg_callbacks.m called with unspecified operation!\n');

return;



function plotmat

global st mancoregvar;

mat = getmat(mancoregvar);

premul=spm_matrix(mat);

set(mancoregvar.hmat_1_1,'string',sprintf('%2.4g',(premul(1,1)) ));
set(mancoregvar.hmat_1_2,'string',sprintf('%2.4g',(premul(1,2)) ));
set(mancoregvar.hmat_1_3,'string',sprintf('%2.4g',(premul(1,3)) ));
set(mancoregvar.hmat_1_4,'string',sprintf('%2.4g',(premul(1,4)) ));
set(mancoregvar.hmat_2_1,'string',sprintf('%2.4g',(premul(2,1)) ));
set(mancoregvar.hmat_2_2,'string',sprintf('%2.4g',(premul(2,2)) ));
set(mancoregvar.hmat_2_3,'string',sprintf('%2.4g',(premul(2,3)) ));
set(mancoregvar.hmat_2_4,'string',sprintf('%2.4g',(premul(2,4)) ));
set(mancoregvar.hmat_3_1,'string',sprintf('%2.4g',(premul(3,1)) ));
set(mancoregvar.hmat_3_2,'string',sprintf('%2.4g',(premul(3,2)) ));
set(mancoregvar.hmat_3_3,'string',sprintf('%2.4g',(premul(3,3)) ));
set(mancoregvar.hmat_3_4,'string',sprintf('%2.4g',(premul(3,4)) ));
set(mancoregvar.hmat_4_1,'string',sprintf('%2.4g',(premul(4,1)) ));
set(mancoregvar.hmat_4_2,'string',sprintf('%2.4g',(premul(4,2)) ));
set(mancoregvar.hmat_4_3,'string',sprintf('%2.4g',(premul(4,3)) ));
set(mancoregvar.hmat_4_4,'string',sprintf('%2.4g',(premul(4,4)) ));


return;


function [VG,VF] = getimg(mancoregvar)

flags.sep = 1;

VG = spm_vol(mancoregvar.targetimage);
if ~isfield(VG, 'uint8')
    VG.uint8 = loaduint8(VG);
    vxg      = sqrt(sum(VG.mat(1:3,1:3).^2));
    fwhmg    = sqrt(max([1 1 1]*flags.sep(end)^2 - vxg.^2, [0 0 0]))./vxg;
    VG       = smooth_uint8(VG,fwhmg);
end

VF = spm_vol(mancoregvar.sourceimage);
if ~isfield(VF, 'uint8')
    VF.uint8 = loaduint8(VF);
    vxg      = sqrt(sum(VF.mat(1:3,1:3).^2));
    fwhmg    = sqrt(max([1 1 1]*flags.sep(end)^2 - vxg.^2, [0 0 0]))./vxg;
    VF       = smooth_uint8(VF,fwhmg);
end


function mat = getmat(mancoregvar)

angl_pitch=get(mancoregvar.hpitch,'Value');
angl_roll=get(mancoregvar.hroll,'Value');
angl_yaw=get(mancoregvar.hyaw,'Value');
dist_x=get(mancoregvar.hx,'Value');
dist_y=get(mancoregvar.hy,'Value');
dist_z=get(mancoregvar.hz,'Value');

mat=[dist_x dist_y dist_z angl_pitch angl_roll angl_yaw 1 1 1 0 0 0];


% ======== Copied from spm_coreg.m ======== %

%==========================================================================
% function udat = loaduint8(V)
%==========================================================================
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 && V.pinfo(1) == 2
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3)
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        img = img(isfinite(img));
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
    end
end

% Another pass to find a maximum that allows a few hot-spots in the data.
nh = 2048;
h  = zeros(nh,1);
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    img = img(isfinite(img));
    img = round((img+((mx-mn)/(nh-1)-mn))*((nh-1)/(mx-mn)));
    h   = h + accumarray(img,1,[nh 1]);
end
tmp = [find(cumsum(h)/sum(h)>0.9999); nh];
mx  = (mn*nh-mx+tmp(1)*(mx-mn))/(nh-1);

% Load data from file indicated by V into an array of unsigned bytes.
udat = zeros(V.dim,'uint8');
st = rand('state'); % st = rng;
rand('state',100); % rng(100,'v5uniform'); % rng('defaults');
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0
        udat(:,:,p) = uint8(max(min(round((img-mn)*(255/(mx-mn))),255),0));
    else
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(max(min(round((img+r-mn)*(255/(mx-mn))),255),0));
    end
end
rand('state',st); % rng(st);


%==========================================================================
% function acc = paccuracy(V,p)
%==========================================================================
function acc = paccuracy(V,p)
if ~spm_type(V.dt(1),'intt')
    acc = 0;
else
    if size(V.pinfo,2)==1
        acc = abs(V.pinfo(1,1));
    else
        acc = abs(V.pinfo(1,p));
    end
end


%==========================================================================
% function V = smooth_uint8(V,fwhm)
%==========================================================================
function V = smooth_uint8(V,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);



%==========================================================================
% function o = optfun(x,VG,VF,s,cf,fwhm)
%==========================================================================
function o = optfun(x,VG,VF,s,cf,fwhm)
% The function that is minimised.
if nargin<6, fwhm = [7 7];   end
if nargin<5, cf   = 'mi';    end
if nargin<4, s    = [1 1 1]; end

% Voxel sizes
vxg = sqrt(sum(VG.mat(1:3,1:3).^2));sg = s./vxg;

% Create the joint histogram
H = spm_hist2(VG.uint8,VF.uint8, VF.mat\spm_matrix(x(:)')*VG.mat ,sg);

% Smooth the histogram
lim  = ceil(2*fwhm);
krn1 = spm_smoothkern(fwhm(1),-lim(1):lim(1)) ; krn1 = krn1/sum(krn1); H = conv2(H,krn1);
krn2 = spm_smoothkern(fwhm(2),-lim(2):lim(2))'; krn2 = krn2/sum(krn2); H = conv2(H,krn2);

% Compute cost function from histogram
H  = H+eps;
sh = sum(H(:));
H  = H/sh;
s1 = sum(H,1);
s2 = sum(H,2);

switch lower(cf)
    case 'mi'
        % Mutual Information:
        H   = H.*log2(H./(s2*s1));
        mi  = sum(H(:));
        o   = -mi;
    case 'ecc'
        % Entropy Correlation Coefficient of:
        % Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
        % "Multimodality image registration by maximisation of mutual
        % information". IEEE Transactions on Medical Imaging 16(2):187-198
        H   = H.*log2(H./(s2*s1));
        mi  = sum(H(:));
        ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
        o   = -ecc;
    case 'nmi'
        % Normalised Mutual Information of:
        % Studholme,  Hill & Hawkes (1998).
        % "A normalized entropy measure of 3-D medical image alignment".
        % in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
        nmi = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
        o   = -nmi;
    case 'ncc'
        % Normalised Cross Correlation
        i     = 1:size(H,1);
        j     = 1:size(H,2);
        m1    = sum(s2.*i');
        m2    = sum(s1.*j);
        sig1  = sqrt(sum(s2.*(i'-m1).^2));
        sig2  = sqrt(sum(s1.*(j -m2).^2));
        [i,j] = ndgrid(i-m1,j-m2);
        ncc   = sum(sum(H.*i.*j))/(sig1*sig2);
        o     = -ncc;
    otherwise
        error('Invalid cost function specified');
end


