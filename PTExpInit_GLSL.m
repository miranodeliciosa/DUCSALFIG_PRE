function [ varargout ] = PTExpInit_GLSL( Scr , onPRPX)
%function [ ScrNum, screensize, xCenter, yCenter, window, framerate  ] = ScreenInit_GLSL( Scr )
%making some standard screen preparations for visual stimulation
%   importantly, GLSL functionality is initialized here, make sure that the
%   graphics card is capable of it, otherwise the function may crash
%   input is a structure with the following fields:
fprintf('PSYCHTOOLBOX: adjusting video settings...requires some seconds')
AssertOpenGL
% PsychDefaultSetup(2); % magic, does: AssertOpenGL, execute KbName('UnifyKeyNames') and Screen('ColorRange', window, 1, [], 1)
% PsychImaging('PrepareConfiguration');% - Prepare setup of imaging pipeline for onscreen window.
% This is the first step in the sequence of configuration steps.

% defaults
if nargin < 2
    onPRPX = 0;
end
if ~isfield(Scr,'PRPXres')
    Scr.PRPXres = [1920 1080]; % this shouldn't be changed as its values can only be set up within propixx windows utility "vputil"
end
if ~isfield(Scr,'RefRate') || onPRPX
    Scr.RefRate = 120; % Hz
end
if ~isfield(Scr,'BckGrCol')
    Scr.BckGrCol = [0 0 0 1]; % black background, fully opaque
end

% get number of available screens
if ~isfield(Scr,'ScrNum')
    Scr.ScrNum = [];
end
screens=Screen('Screens');
if isempty(Scr.ScrNum) || all(screens == 0)
    Scr.ScrNum=max(screens);
end
if ~isfield(Scr,'PRPXmode')
    Scr.PRPXmode = 0;
end
% Propixx
if onPRPX
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', Scr.PRPXmode); % 2 for 480, 5 for 1440 Hz, 0 for normal (120Hz)
    Datapixx('RegWrRd');
end

% Query size of screen: 
% [screenXpixels, screenYpixels] = Screen('WindowSize', window); 
screensize=Screen('Rect', Scr.ScrNum);

% some keyboard settings
KbName('UnifyKeyNames')
keymap = zeros(1,256); % the keymap index vector for KbQueue function family

% adjust screen parameters if necessary, only possible on linux
if IsLinux 
    if ~all(screensize(1,3:4) == Scr.PRPXres)
        Screen('ConfigureDisplay', 'Scanout', Scr.ScrNum, 0, Scr.PRPXres(1),Scr.PRPXres(2),Scr.RefRate); %'0' is the virtual screen output id?? sets up the target screen with prompted parameters
        WaitSecs(10); % allow screen changes to take place
        screensize=Screen('Rect', Scr.ScrNum);
    end
    % parallel port
    try ppdev_mex('Close', 1); % open triggerport
    catch me
        disp(me)
    end
    try ppdev_mex('Open', 1); % open triggerport
    catch me
        disp(me)
    end
    lptwrite(1,0); % ensure it's set to zero (solves the LPT-4-issue!)
%     RespDev = -1;
        RespDev = GetKeyboardIndices('Logitech USB Keyboard',[],3);
        ExpDev = GetKeyboardIndices('Chicony HP Business Slim Keyboard',[],3);
%     RespDev = GetKeyboardIndices('Chicony HP Elite USB Keyboard',[],3); % either enter 'product' (name) | 'serialNumber' | locationID
else
    RespDev = -1; % for KbQueueX functions this will indicate to use the default device
end

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(screensize);

% psychtoolbox start screen in black (instead of white)
Screen('Preference', 'VisualDebugLevel', 1);
% Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SkipSyncTests', 2);
% Screen('Preference', 'SkipSyncTests', 1);
% open a psychtoolbox window in background color
window = Screen('OpenWindow', Scr.ScrNum, Scr.BckGrCol.*255);%, [0,0,800,600]); % 1 = main screen (0|2 -> main screen with menu bar|second monitor), color = 0 (white|[1 1 1 a]), rect = default (fullscreen|[0,0,width,height])
%window = Screen('OpenWindow', Scr.ScrNum, Scr.BckGrCol.*255, [0 0 500 500]);%, [0,0,800,600]); % 1 = main screen (0|2 -> main screen with menu bar|second monitor), color = 0 (white|[1 1 1 a]), rect = default (fullscreen|[0,0,width,height])
% window = PsychImaging('OpenWindow',ScrNum,BckGrCol);
% Screen('ColorRange', window,1); % clamping color range to 1 also has good effects on graphic card communication...sometimes
Screen('ColorRange', window,1,[],1); % clamping color range to 1 also has good effects on graphic card communication...sometimes
% set blending modes
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Make sure the GLSL shading language is supported:
AssertGLSL;

framerate=Screen('NominalFramerate', window); % in Hz
if framerate ~= Scr.RefRate
    warning('measured screen refresh rate of %s does not correspond to queried rate of %s...presentation assumes the latter.',num2str(framerate),num2str(Scr.RefRate))
    framerate= Scr.RefRate;
end

% block button presses to matlab window STRG + C to exit and ListenChar(0) at end
ListenChar(-1) 
%ListenChar(0) 
% set priority to high (1)
Priority(1); 

% structure output
varargout{1} = Scr.ScrNum;
varargout{2} = screensize;
varargout{3} = xCenter;
varargout{4} = yCenter;
varargout{5} = window;
varargout{6} = framerate;
varargout{7} = RespDev;
varargout{8} = keymap;
varargout{9} = Scr.RefRate;

fprintf('done!\n')
end

