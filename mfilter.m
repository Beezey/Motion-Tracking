clear all
close all
clc
imtool close all



folder = fileparts(which('classroom.mp4'));
movieFullFileName = fullfile(folder, 'classroom.mp4');

path = '..\ExtractMovieFrames';
o_f = 'Output';
out_folder = fullfile(path,o_f)

videoObject = VideoReader(movieFullFileName);
% Determine how many frames there are.
numberOfFrames = videoObject.NumberOfFrames;
vidHeight = videoObject.Height;
vidWidth = videoObject.Width;

im_start = 1;
im_end = numberOfFrames;



num_images = im_end-im_start;
for i = im_start:im_end
    thisFrame = read(videoObject, i);
%     images(:,:,i+1) = rgb2gray(imread(sprintf('Frame %4d of %d.', i, numberOfFrames)));
    images(:,:,i+1) = rgb2gray(thisFrame);
end

filter = [-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'Simple_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
%simple 1D
% 
filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'Simple_1D_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
% gaussian
tsigma = 1.2;
gfilter = fspecial('gaussian',[1,5],tsigma);
filter = gradient(gfilter);
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 5;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'Gaussian_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
% higher tsigma

tsigma = 1;
gfilter = fspecial('gaussian',[1,round(5*tsigma)],tsigma);
filter = gradient(gfilter);
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),...
    num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 5;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'higher_tsigma_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
%3x3 box filter
box_size = 3;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i),...
        fspecial('average',[box_size,box_size]));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1),...
    size(images_box(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_box(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_box;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'3x3_box_Spatial_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
%5x5 box filter
box_size = 5;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i),...
        fspecial('average',[box_size,box_size]));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1),...
    size(images_box(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_box(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_box;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'5x5_box_Spatial_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
%2D Gaussian
tsigma = 1;
for i = 1:num_images
    images_gaussian(:,:,i) = imfilter(images(:,:,i),...
        fspecial('gaussian',[5,5],tsigma));
end

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images_gaussian(:,:,1),1),...
    size(images_gaussian(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images_gaussian(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images_gaussian(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end
threshold = 15;
mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images_gaussian;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'2D_Gaussian_filter.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)
%%
%threshold optimization

filter = 0.5*[-1 0 1];
filter_depth = (length(filter)-1)/2;
filtered_image = zeros(size(images(:,:,1),1),...
    size(images(:,:,1),2),num_images);
mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    filtered_values = zeros(size(images(:,:,1)));
    for l = 1:length(filter)
        filtered_values = filtered_values + ...
            double(images(:,:,i-filter_depth+l-1))*filter(l);
    end
    filtered_image(:,:,i) = abs(filtered_values);
end

sub_images = double(images(size(images,1)-20:size(images,1),1:20,:));
E = 1/num_images * sum(sub_images,3);
% E mean squared difference
E_msd = 0;
for i = 1:num_images
    E_msd = E_msd + (E - sub_images(:,:,i)).^2;
end
sigma = sqrt(1/(num_images - 1)*E_msd);
est_sigma = mean(mean(sigma));
threshold = 3*est_sigma;

mask(filtered_image >= threshold) = 1;
mask(filtered_image < threshold) = 0;

result = images;
result(mask==1) = 255;

implay(result)

outputVideo = VideoWriter(fullfile(out_folder,'Threshold_Optimization.avi'));
outputVideo.FrameRate = videoObject.FrameRate;
open(outputVideo)

for ii = 1:length(result)
   img = result(:,:,ii);
   writeVideo(outputVideo,img)
end

close(outputVideo)