function spm_reorient(ImagesFiles2Process, M)
%spm_reorient(ImagesFiles2Process, M)
% Just applies a transformation matrix to a bunch of images (only changes the header)
% ImagesFiles2Process{n,1} : a cell of n images fullpath filenames
% M - reorientation matrix that will move the source image to the target
% image

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

spm_get_space(ImagesFiles2Process{1}, M*spm_get_space(ImagesFiles2Process{1}));

for FileInd=1:size(ImagesFiles2Process,1)
    spm_progress_bar('Set',FileInd);
    spm_get_space(ImagesFiles2Process{FileInd,1}, M*Mats(:,:,FileInd));
end

spm_progress_bar('Clear')

end

