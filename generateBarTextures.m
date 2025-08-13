function [barTex, dstRects, angles] = generateBarTextures(ps, SBA, SBAin)
% generateBarTextures - Prepares bar stimulus textures, positions, and rotation angles for Psychtoolbox
%
% This function generates rectangular bar textures for use with Screen('DrawTextures') in Psychtoolbox.
% The bars are evenly distributed in a grid within a fixed area centered on the screen.
% Multiple shape events per trial are supported; each can trigger rotation of a subset of bars at specific frames.
%
% USAGE:
%   [barTex, dstRects, angles] = generateBarTextures(ps, SBA, SBAin)
%
% INPUTS:
%   ps               - Structure containing Psychtoolbox window handle (ps.window) and framerate (ps.framerate)
%   SBA              - Structure containing stimulus parameters:
%                      - numBars: [cols, rows]
%                      - sizeBars: [length, width]
%                      - defaultAngle: default bar orientation
%                      - anglesToRotate: vector of rotation angles
%                      - colors: n x 3 RGB matrix (0–1 values)
%                      - freq: color-switching frequency
%   SBAin            - Trial-specific structure with:
%                      - trial.frames: total number of frames
%                      - trial.event.shape: vector of shape indices (which bar groups to rotate)
%                      - trial.event.onset: vector of onset frames (same length as shape)
%                      - trial.event.duration: scalar (duration of each event in sec)
%                      - trial.event.contrast: vector of contrast (index into anglesToRotate)
%
% OUTPUTS:
%   barTex           - Cell array of bar textures per frame
%   dstRects         - Destination rectangles (4 × totalBars)
%   angles           - Matrix (totalBars × numFrames) with rotation angle (in degrees) per bar and frame

% screen parameters
window = ps.window;
frameRate = ps.framerate;

% stimulus parameters
numBars = SBA.numBars;
sizeBars = SBA.sizeBars;
defaultAngle = SBA.defaultAngle;
anglesToRotate = SBA.event.anglesToRotate;
directionToRotate = SBAin.trial.event.direction;
colors = SBA.colors;
freq = SBA.freq;
numCol = size(SBA.colors,1);

% event parameters
numFrames = SBAin.trial.frames;
onset = SBAin.trial.event.onset;
if any(isnan(onset)); onset(isnan(onset)) = 1; end 
contrast = SBAin.trial.event.contrast; 
if any(isnan(contrast)); contrast(isnan(contrast)) = 1; end 
shapeVec = SBAin.trial.event.shape; 
if any(isnan(shapeVec)); shapeVec(isnan(shapeVec)) = 1; end 

% define barsToRotate for all possible shapes
[iCols, iRows] = extractColorIndices('images', [0, 0, 0], SBA.numBars);  % returns index maps
allBarsToRotate = cell(numel(iCols), 1);
for s = 1:numel(iCols)
    allBarsToRotate{s} = sub2ind([SBA.numBars(1), SBA.numBars(2)], iRows{s}, iCols{s});
end

% define color switching over frames
vec = []; 
dutycycle = 0.5;
for i=1:ceil(numFrames/(numCol*dutycycle*frameRate/freq))
    vec = [vec randperm(numCol)];
    while sum(diff(vec)==0)>0
        vec(end-(numCol-1):end) = randperm(numCol); 
    end
end
ColChange = reshape(repmat(vec, 0.5*frameRate/freq,1), [], 1)';

% constants
areaSize = SBA.size(1);
barLength = sizeBars(1);
barWidth  = sizeBars(2);
nCols = numBars(1);
nRows = numBars(2);
totalBars = nCols * nRows;

% compute grid
spacingX = areaSize / (nCols + 1);
spacingY = areaSize / (nRows + 1);
[centerX, centerY] = RectCenter(Screen('Rect', window));
gridWidth  = (nCols - 1) * spacingX;
gridHeight = (nRows - 1) * spacingY;
xOffset = centerX - gridWidth / 2;
yOffset = centerY - gridHeight / 2;

positions = zeros(totalBars, 2);
barIndex = 1;
for row = 1:nRows
    for col = 1:nCols
        x = (col - 1) * spacingX + xOffset;
        y = (row - 1) * spacingY + yOffset;
        positions(barIndex, :) = [x, y];
        barIndex = barIndex + 1;
    end
end

% generate textures
barTextures = cell(1, numCol);
for c = 1:numCol
    color = uint8(reshape(colors(c,:)*255, 1, 1, 3));
    barImage = repmat(color, barWidth, barLength, 1);
    barTextures{c} = Screen('MakeTexture', window, barImage);
end

% assign textures to frames
barTex = cell(1, numFrames);
for iFrame = 1:numFrames
    barTex{iFrame} = barTextures{ColChange(iFrame)};
end 

% destination rectangles
dstRects = zeros(4, totalBars);
for i = 1:totalBars
    cx = positions(i, 1);
    cy = positions(i, 2);
    rect = CenterRectOnPointd([0 0 barLength barWidth], cx, cy);
    dstRects(:, i) = rect';
end

% initialize angle matrix
angles = defaultAngle * ones(totalBars, numFrames);

% apply event-wise rotations
eventDurFrames = SBAin.trial.event.duration * frameRate;

for i_ev = 1:length(shapeVec)
    shapeIdx = shapeVec(i_ev);
    barIndices = allBarsToRotate{shapeIdx};  % linear indices

    % jitter exact shape position
    sz = [SBA.numBars(1), SBA.numBars(2)];
    [r, c] = ind2sub(sz, barIndices);
    
    dr = randi([-1, 1], 1);
    dc = randi([-1, 1], 1);
    
    wrap = @(x,L) mod(x-1, L) + 1;  % 1-based wrap
    r_jit = wrap(r + dr, sz(1));
    c_jit = wrap(c + dc, sz(2));
    
    jitteredBarIndices = sub2ind(sz, r_jit, c_jit);
    targetBarIndices = jitteredBarIndices;
    
    % compute rotation offset
    if contrast(i_ev) > 0
        rotationOffset = anglesToRotate(contrast(i_ev))*directionToRotate(i_ev);
    else
        rotationOffset = 0;
    end

    % apply rotation offset to jittered barsIndices
    evFrames = onset(i_ev):(onset(i_ev) + eventDurFrames - 1);
    evFrames = evFrames(evFrames <= numFrames); % prevent overflow
    angles(targetBarIndices, evFrames) = angles(targetBarIndices, evFrames) + rotationOffset;
end

end


% %function [barTextures, dstRects, angles] = generateBarTextures(window, numBars, barsToRotate, sizeBars, defaultAngle, anglesToRotate, colors)
% function [barTex, dstRects, angles] = generateBarTextures(ps, SBA, SBAin)
% % generateBarTextures - Prepares bar stimulus textures, positions, and rotation angles for Psychtoolbox
% %
% % This function generates rectangular bar textures for use with Screen('DrawTextures') in Psychtoolbox.
% % The bars are evenly distributed in a grid within a *areasize* × *areasize* pixel area and centered on the display window.
% % A subset of bars can be rotated by specified angles relative to a default orientation.
% %
% % USAGE:
% %   [barTextures, dstRects, angles] = generateBarTextures(window, numBars, numBarsToRotate, ...
% %       barsToRotate, sizeBars, defaultAngle, anglesToRotate, colors)
% %
% % INPUTS:
% %   window           - Psychtoolbox window pointer from Screen('OpenWindow').
% %   numBars          - [cols, rows]; defines the grid layout of bars, cols x rows = total number of bars
% %   numBarsToRotate  - Number or proportion of bars to rotate. Ignored if barsToRotate is provided.
% %   barsToRotate     - (Optional) cell array of n×2 matrix of [col, row] positions (MATLAB-style) of bars to rotate.
% %                      If empty or not provided, bars are selected randomly based on numBarsToRotate.
% %   sizeBars         - [length, width] of each bar in pixels (pre-rotation).
% %   defaultAngle     - Default orientation of all bars in degrees (e.g., 0 = horizontal).
% %   anglesToRotate   - Vector of angles (in degrees) to be added to selected bars' default orientation.
% %                      One output set (angles) is created per angle in this vector.
% %   colors           - Cell array of RGB vectors (e.g., {[255 0 0]; [0 255 0]}); one texture per color.
% %
% % OUTPUTS:
% %   barTextures      - Cell array of texture pointers (length = number of colors).
% %   dstRects         - 4×n matrix of destination rectangles for each bar.
% %   angles           - nxnflips matrix of rotation angles in degrees, where n = total number of bars, and nflips = number of
% %                      flips
% 
% %
% % NOTE:
% %   This function assumes a fixed grid area of *areasize* × *areasize* px centered in the window. All bars are
% %   positioned based on their grid index and the area is centered on screen.
% %
% % EXAMPLE USAGE:
% %   Screen('DrawTextures', window, barTextures{1}, [], dstRects, angles(2, :));
% %   Screen('Flip', window);
% 
% % screen parameters
% window = ps.window;
% frameRate = ps.framerate;
% 
% 
% % stimulus parameters
% numBars = SBA.numBars;
% sizeBars = SBA.sizeBars;
% defaultAngle = SBA.defaultAngle;
% anglesToRotate = SBA.event.anglesToRotate;
% colors = SBA.colors;
% freq = SBA.freq;
% numCol = size(SBA.colors,1);
% 
% % event parameters
% numFrames = SBAin.trial.frames;
% onset = SBAin.trial.event.onset; if isnan(onset); onset=1;end 
% contrast = SBAin.trial.event.contrast; if isnan(contrast); contrast=1;end 
% shape = SBAin.trial.event.shape; if isnan(shape); shape=1;end % index for barstToRotate
% 
% % define bars to rotate
% [iCols, iRows] = extractColorIndices('images', [0, 0, 0], SBA.numBars);
% barsToRotate = {iRows{shape}; iCols{shape}}; 
% 
% % define frame structure
% vec = []; 
% dutycycle = 0.5;
% for i=1:ceil(numFrames/(numCol*dutycycle*frameRate/freq))
%     vec = [vec randperm(numCol)]; 
%     %vec = [vec 1:numCol]; 
%     while sum(diff(vec)==0)>0
%         vec(end-(numCol-1):end) = randperm(numCol); 
%     end
% end
% ColChange = reshape(repmat(vec, 0.5*frameRate/freq,1), [], 1)'; % GENERATEBARTEXTURE
% %figure; plot(ColChange)
% 
% % Constants
% areaSize = 500;  % size of the grid area in pixels
% barLength = sizeBars(1);
% barWidth  = sizeBars(2);
% nCols = numBars(1);
% nRows = numBars(2);
% totalBars = nCols * nRows;
% 
% % Compute grid spacing
% spacingX = areaSize / (nCols + 1);
% spacingY = areaSize / (nRows + 1);
% 
% % Compute center of the window
% windowRect = Screen('Rect', window);
% [centerX, centerY] = RectCenter(windowRect);
% 
% % Compute total grid size
% gridWidth  = (nCols - 1) * spacingX;
% gridHeight = (nRows - 1) * spacingY;
% 
% % Compute top-left origin so grid is centered
% xOffset = centerX - gridWidth / 2;
% yOffset = centerY - gridHeight / 2;
% 
% % Generate all bar center positions
% positions = zeros(totalBars, 2); % [x, y] per row
% barIndex = 1;
% for row = 1:nRows
%     for col = 1:nCols
%         x = (col - 1) * spacingX + xOffset;
%         y = (row - 1) * spacingY + yOffset;
%         positions(barIndex, :) = [x, y];
%         barIndex = barIndex + 1;
%     end
% end
% 
% % Generate base texture(s)
% numColors = size(colors,1);
% barTextures = cell(1, numColors);
% 
% for c = 1:numColors
%     color = uint8(reshape(colors(c,:)*255, 1, 1, 3));
%     barImage = repmat(color, barWidth, barLength, 1);
%     barTextures{c} = Screen('MakeTexture', window, barImage);
% end
% % assign bartextures to frames
% barTex = cell(1, numFrames);
% for iFrame = 1:numFrames
%     barTex{iFrame} = barTextures{ColChange(iFrame)};
% end 
% 
% % Generate destination rectangles (all same size, centered on grid)
% dstRects = zeros(4, totalBars);
% for i = 1:totalBars
%     cx = positions(i, 1);
%     cy = positions(i, 2);
%     rect = CenterRectOnPointd([0 0 barLength barWidth], cx, cy);
%     dstRects(:, i) = rect';
% end
% 
% % Determine rotated bars
% if exist('barsToRotate', 'var') && ~isempty(barsToRotate)
%     % Convert [col, row] pairs to linear indices
%     linearIdx = sub2ind([nCols, nRows], barsToRotate{1}, barsToRotate{2});
% else
%     if numBarsToRotate < 1
%         numToRotate = round(numBarsToRotate * totalBars);
%     else
%         numToRotate = round(numBarsToRotate);
%     end
%     linearIdx = randperm(totalBars, numToRotate);
% end
% 
% % Generate angle matrix
% nAngles = numel(anglesToRotate);
% %angles = defaultAngle*ones(nAngles, totalBars, numFrames);
% angles = defaultAngle*ones(totalBars, numFrames);
% eventFrames = [];
% for i_ev = 1:length(SBAin.trial.event.onset)
%     eventFrames = SBAin.trial.event.onset(i_ev):(SBAin.trial.event.onset(i_ev) + SBAin.trial.event.duration*frameRate);
%     for iFrame = eventFrames
%         if SBAin.trial.event.contrast>0
%         rotationOffset = anglesToRotate(SBAin.trial.event.contrast(i_ev));
%         else 
%             rotationOffset=0;
%         end 
%         angles(linearIdx,iFrame) = angles(linearIdx,iFrame) + rotationOffset;
%     end 
%     eventFrames = [];
% end 
% 
% % 
% % for iAngle = 1:nAngles
% %     rotationOffset = anglesToRotate(iAngle);
% %     angles(iAngle, :) = defaultAngle;  % initialize all bars to default
% %     angles(iAngle, linearIdx) = defaultAngle + rotationOffset;
% % end
% 
% end
