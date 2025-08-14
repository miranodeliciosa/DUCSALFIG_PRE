function [iRows, iCols] = extractColorIndices(imagesFolder, color, sizeToReshape)
% extractColorIndices - Load and process images to extract pixel indices of a specific RGB color
%
% INPUTS:
%   imagesFolder   - Name of the folder containing PNG images (relative to current folder)
%   color          - 1x3 RGB vector specifying the target color (e.g., [255 0 0])
%   sizeToReshape  - 1x2 vector [height, width] for resizing each image
%
% OUTPUTS:
%   iRows          - Cell array; each cell contains row indices of matching pixels in one image
%   iCols          - Cell array; each cell contains column indices of matching pixels in one image
%
% NOTES:
%   - Alpha channels are ignored
%   - RGB values must match exactly
%   - Images are assumed to be RGB and resized using imresize

% Get all PNG files in the folder
%imageFiles = dir(fullfile(imagesFolder, '*.png'));
imageFiles = dir(fullfile(imagesFolder, '*.tiff'));

% Initialize output
iRows = cell(1, numel(imageFiles));
iCols = cell(1, numel(imageFiles));

% Process each image
for k = 1:numel(imageFiles)
    % Read image
    imgPath = fullfile(imagesFolder, imageFiles(k).name);
    img = imread(imgPath);
    % debugging 
    
    % Ignore alpha channel if present
    if size(img, 3) > 3
        img = img(:, :, 1:3);
    end

    % Resize image
    imgResized = imresize(img, sizeToReshape);

    % Find pixels matching the specified color
%     if isempty(img(img(:)~=0 & img(:)~=255))
%         mask = imgResized(:, :, 1) == color(1) & ...
%                imgResized(:, :, 2) == color(2) & ...
%                imgResized(:, :, 3) == color(3);
%     else
       mask =  imgResized(:, :, 1) <= 180 & ...
                imgResized(:, :, 2) <= 180 & ...
                imgResized(:, :, 3) <= 180;
      %unique(img(img(:)~=0 & img(:)~=255))
%     end 

    % Extract row and column indices
    [rows, cols] = find(mask);
    iRows{k} = rows;
    iCols{k} = cols;
end
end
