function FixTex = MakeFixationTextures(ps, p)
    % Create a fixation cross texture centered on the screen

    window = ps.window;
    center = [ps.xCenter ps.yCenter];  % should be [x y]
    color = p.crs.color;
    imgSize = p.crs.size;
    lineWidth = p.crs.width;

    % Cross dimensions (both arms fit in imgSize)
    texSize = imgSize;  % size of the square image

    % Horizontal and vertical bar definitions, centered at 0,0
    hLine = [-imgSize/2, -lineWidth/2, imgSize/2, lineWidth/2];
    vLine = [-lineWidth/2, -imgSize/2, lineWidth/2, imgSize/2];

    % Offscreen window (square)
    crossMat = Screen('OpenOffscreenWindow', window, 255*[.05 .05 .05], [0 0 texSize texSize]);

    % Center coordinates of offscreen texture
    [cx, cy] = RectCenter([0 0 texSize texSize]);

    % Center the bars inside the offscreen window
    hLine = CenterRectOnPoint(hLine, cx, cy);
    vLine = CenterRectOnPoint(vLine, cx, cy);

    % Draw bars
    Screen('FillRect', crossMat, ps.input.BckGrCol, [0 0 imgSize imgSize]); % rect = [left top right bottom];
    Screen('FillRect', crossMat, color, hLine);
    Screen('FillRect', crossMat, color, vLine);

    % Convert to image matrix, then to texture
    imageMatrix = Screen('GetImage', crossMat);
    FixTex = Screen('MakeTexture', window, imageMatrix);

    % Optional test: draw it
%     imagesc(imageMatrix)
    % screenRect = Screen('Rect', window);
    % dstRect = CenterRectOnPoint([0 0 texSize texSize], center(1), center(2));
    % Screen('DrawTexture', window, FixTex, [], dstRect);
    % Screen('Flip', window);


% TODO: add black box that equals the size of the bars to the fixation box
% ODER: lasse den zentralen weg... 

end
