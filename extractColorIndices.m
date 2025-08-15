function [iRows, iCols] = extractColorIndices(imagesFolder, color, sizeToReshape)
% extractColorIndices - Extracts pixel locations of shape masks from image files
%
%   [iRows, iCols] = extractColorIndices(imagesFolder, color, sizeToReshape)
%
%   This function reads shape mask images from a given folder and identifies
%   the pixel locations that match a specified RGB color or fall below a luminance
%   threshold. The resulting pixel positions define bar locations for shape events
%   in the Tiltanic experiment.
%
%   INPUTS:
%       imagesFolder   - Path to folder containing shape images (*.tiff format)
%                        Images should be in RGB format. Alpha channels are ignored.
%       color          - RGB color vector [R G B] (0–255) indicating which color to detect
%                        Currently not strictly applied — instead, a low luminance threshold
%                        is used to detect figure regions.
%       sizeToReshape  - Desired [height, width] in pixels for resizing each image,
%                        typically matching the number of SBA bars (e.g., [12, 12])
%
%   OUTPUTS:
%       iRows          - Cell array (1 × nShapes); row indices of figure pixels per image
%       iCols          - Cell array (1 × nShapes); column indices of figure pixels per image
%
%   FUNCTION BEHAVIOR:
%       • Each shape image is resized to the SBA grid size
%       • Black/dark regions (RGB ≤ 180) are interpreted as figure shapes
%       • The positions of dark pixels are extracted to determine which SBA bars
%         will rotate during shape events in `generateBarTextures`
%
%   ASSUMPTIONS & NOTES:
%       • Intended for use with grayscale or binary mask-like TIFF images
%       • RGB matching is currently implemented via a luminance threshold, not exact match
%       • Output format aligns with SBA grid indexing (rows/columns for bar placement)
%
%   DEPENDENCIES:
%       Used by: generateBarTextures.m
%
%   SEE ALSO:
%       generateBarTextures, run_Tiltanic
%
%   Author: Sebastian Wehle, Leipzig (2025)


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
