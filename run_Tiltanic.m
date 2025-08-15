function [] = run_Tiltanic(sub,targetshape, flag_training)
% run_Tiltanic - Runs the SSVEP-based 2AFC figure-ground segmentation experiment
%
%   run_Tiltanic(sub, targetshape, flag_training)
%
%   This script presents a texture segmentation task designed to estimate
%   the perceptual threshold for shape detection based on orientation contrast.
%   Participants view a superimposed RDK and static background bar array (SBA),
%   where brief events rotate a subset of bars to form geometric shapes.
%   The participant’s task is to indicate whether the perceived shape matches
%   a predefined target shape in a two-alternative forced-choice (2AFC) task.
%
%   INPUTS:
%       sub           - Subject ID (numeric). Used for log file naming and random seed
%       targetshape   - Index of the target shape (integer from 1–5)
%                          1 = 'mask'
%                          2 = 'dirttruck'
%                          3 = 'hat'
%                          4 = 'sumo'
%                          5 = 'cactus'
%       flag_training - 1 to run training with shape preview and practice trials
%                     - 0 to skip training and run the experiment directly
%
%   USAGE EXAMPLE:
%       run_Tiltanic(99, 3, 1);    % Runs training and experiment for subject 99 with 'hat' as target shape
%
%   EXPERIMENT STRUCTURE:
%       • Training phase (optional):
%           - Shape familiarization: each shape is shown once with text prompts
%           - Practice trials with reduced difficulty (fewer events)
%           - Repetition allowed until participant confirms readiness
%       • Experiment phase:
%           - 16 trials (1.5–2 min each) with 40 randomized segmentation events per trial
%           - Participants respond whether the shape matches the target (2AFC)
%           - Feedback on valid responses is given after each trial
%
%   STIMULI OVERVIEW:
%       • SBA (background bar array):
%           - 12x12 grid of colored bars, flickering at 7.5 Hz
%           - Each trial contains brief segmentation events: subsets of bars are rotated
%             to form distinct shapes (e.g., diamond, hat, sumo, etc.)
%           - Rotation angles correspond to contrast conditions used for threshold estimation
%           - The spatial location and identity of shapes vary across events
%           - Textures for all frames are precomputed using `generateBarTextures`, which:
%               * Generates Psychtoolbox textures for each frame
%               * Computes the spatial positions and rotation angles of each bar
%               * Encodes shape identity, rotation angle (contrast), and onset timing
%
%       • RDK (foreground motion):
%           - Superimposed random dot kinematogram with grey squares
%           - Dots move in random or coherent directions (task-irrelevant motion events)
%           - Motion events are time-jittered and temporally separated from SBA events
%
%   DEPENDENCIES:
%       • generateBarTextures.m: Generates frame-by-frame SBA textures and rotation angles per trial
%       • rand_Tiltanic.m: Randomizes condition assignments and event onsets for SBA and RDK using generate_event_onset_continuous
%       • pres_Tiltanic.m: Runs the full stimulus presentation and response recording for each trial
%       • pres_feedback_Tiltanic.m: Provides feedback after each trial
%       • PTExpInit_GLSL.m: Screen initialization and hardware configuration

%
%   RESPONSES:
%       • 'j' key - target shape was shown
%       • 'f' key - non-target shape was shown
%       • Responses valid only if given within 200–1000 ms after shape onset
%
%   DATA LOGGING:
%       • Saves trial-wise button presses, timings, and stimulus parameters
%       • Logfiles are saved per subject to platform-specific paths
%
%   AUTHORS:
%       Originally adapted from FShiftBase by M. Dotzer & C. Gundlach (2023)
%       Modified and extended by Sebastian Wehle (2025, Leipzig)
%
%   CHANGES IN THIS VERSION:
%       • Static background bar arrays with event-based rotation
%       • Integrated single RDK stream (foreground motion, task-irrelevant)
%       • 2AFC task without central fixation task
%       • continous long trial stimulation
%       • get rid of block structure
%       • no isoluminance adjustment
%       • adjusted randomization and response acquisition
%       • instructive shape presentation 
%

if nargin < 3
    help run_Tiltanic
    return
end

%% parameters
% sub = 97; flag_training = 1; flag_block = 1; 

% design 
p.sub                   = sub;                  % subject number
p.flag_training         = flag_training;        % do training

% screen
p.scr_num               = 1;                    % screen number
if IsOSX; p.scr_num       = 0; end              % screen number on mac

p.scr_res               = [1920 1080];          % resolution
p.scr_refrate           = 120;                  % refresh rate in Hz (e.g. 85)
p.scr_color             = [0.05 0.05 0.05 1];   % default: [0.05 0.05 0.05 1]; ; color of screen [R G B Alpha]
p.scr_imgmultipl        = 1;                    %1: 120 Hz, 4: 480Hz

% stimplan %
p.stim.contrast        = [1 2 3 4 5 6 7 8];  % contrasts (1=contrast 1 (=0), 2=contrast 2, 3= contrast 3, 4 = contrast 4 (highest))
p.stim.condition       = [1 2];              % 1 = target shape, 2=non-target shape
p.stim.eventsPerTrial_e  = 40;               % Events per Trial - Experiment
p.stim.eventsPerTrial_t  = 20;               % Events per Trial - Training
p.stim.con_repeats_e     = 40;               % Event per Contrast Condition (SBA) nEventsPerCondition - Experiment
p.stim.con_repeats_t     = 20;               % Event per Contrast Condition (SBA) nEventsPerCondition - Training
    
% text 
p.stim.text             = {'Pause';'Drücken Sie die Leertaste, um fortzufahren.'};

%IDEA: change pre-cue event to background Event, SBA_event
% SBA event
p.stim.precue_event.num         = [0 0 0 0 1 2 3 4];   % ratio of no SBA-events (0) and SBA-events(1,2,3,4), 1:4  contrast; 
p.stim.precue_event.targets     = [1 2 3 4];%[1 3];     % all capture events are distractor events
p.stim.precue_event.length      = 0.125;      % length of precue-event in s              
p.stim.precue_event.min_onset   = 2.25;       % min time before precue-event onset in s
p.stim.precue_event.min_offset  = 2.25;       % min offset from precue-event end to end of trial in s

% RDK and SBA events
p.stim.event.type           = 1;              % types of events (1 = targets only, 2 = targets + distrators)
p.stim.event.duration         = [0.3 1/7.5];    % duration of [RDK; SBA events in s]
p.stim.event.SOA_ms         = [1 4];          % min distance between events in sec
p.stim.event.min_offset     = [1 1];          % min offset from event end to end of trial in s
p.stim.event.ratio          = [1/2];          % ratio of fore- and background event number
p.stim.trialTime_e          = mean(p.stim.event.SOA_ms)*p.stim.eventsPerTrial_e;              % Preset trialtime --> can be adjusted in generate_event_onset_continuous 
p.stim.trialTime_t          = mean(p.stim.event.SOA_ms)*p.stim.eventsPerTrial_t;              % Preset trialtime --> can be adjusted in generate_event_onset_continuous

p.stim.frames_trialTime_t     = p.stim.trialTime_t*p.scr_refrate;
p.stim.frames_trialTime_e     = p.stim.trialTime_e*p.scr_refrate;


p.targ_respwin          = [200 1000];       % time window for responses in ms
% 
% p.colors                = [0.85 .85 .85 1; % define RDK color (light grey 0.7407 0.7407 0.7407)
%                             0    0.8000    0.6980 1;% Bar Col 1
%                             1.0000    0.3137         0 1;% Bar Col 2
%                             0.2745    0.4118    1.0000 1;% Bar Col 3
%                             0.7843    0.3137    1.0000 1;% Bar Col 4
%                             0.6980    0.9020         0 1]; % Bar Col 5]

p.colors                = [0.85 .85 .85 1;                  % RDK color
                        0    0.64    0.5584 1;              % Bar Col 1
                            0.8    0.25096  0 1;            % Bar Col 2
                            0.2196    0.32944    0.80 1;    % Bar Col 3
                            0.62744    0.25096    0.8 1;    % Bar Col 4
                            0.5584    0.7216    0 1];       % Bar Col 5

                        
% RDK FOREGROUND
RDK.RDK.size             = [750 750];                % width and height of RDK in pixel; only even values [720 = 11.2°; 340 = 10.4°] 
RDK.RDK.centershift      = [0 0];                    % x and y deviation from center in pixel
RDK.RDK.col              = [p.colors(1,:); p.scr_color(1:3) 0];% "on" and "off" color
RDK.RDK.freq             = p.scr_refrate/7;          % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK.mov_freq         = 120;                      % Defines how frequently the dot position is updated; 0 will adjust the update-frequency to your flicker frequency (i.e. dot position will be updated with every "on"-and every "off"-frame); 120 will update the position for every frame for 120Hz or for every 1. quadrant for 480Hz 
RDK.RDK.num              = 150;                      % number of dots 85 120
RDK.RDK.mov_speed        = 1;                        % movement speed in pixel
RDK.RDK.mov_dir          = [0 1; 0 -1; -1 0; 1 0];   % movement direction  [0 1; 0 -1; -1 0; 1 0] = up, down, left, right; !!hard coded with 4 directions in rand-scritp
RDK.RDK.dot_size         = 14;                       % 12 = 0.38°; 10 = 0.31°; 8 = 0.25° --- R120 Values!
RDK.RDK.shape            = 1;                        % 1 = square RDK; 0 = ellipse/circle RDK;

RDK.event.type              = 'globalmotion';               % event type global motion
RDK.event.duration          = p.stim.event.duration(1);     % duration of coherent motion event
RDK.event.coherence         = .4;                           % percentage of coherently moving dots
RDK.event.direction         = RDK.RDK(1).mov_dir;           % movement directions for events

% Static background bar array SBA
SBA.size                 = RDK.RDK.size+50;                     % size of bar grid ; +50 looks better on edges with RDK dots
SBA.colors               = p.colors(2:end,1:3);                 % colors of bars
SBA.numBars              = [12 12];                             % prod(SBA.numbars) total; chose even number for uncovered fixation cross
SBA.sizeBars             = [40, 15];                            % bar size ein pxl, if too large bars overlap!
SBA.freq                 = 7.5;                                 % freq of color change on SBA | refrate/freq should be even 
SBA.defaultAngle         = 90;                                  % 90 = upright and 0 = horizontal if size [max min]; 
SBA.event.anglesToRotate_t       = [32:2:46];                   % Rotation angles in training  

SBA.event.anglesToRotate_e       = [10 16 18 20 22 24 35 45];   % around 20°; Rotation angles in experiment  
%SBA.event.anglesToRotate_e       = [1 2 3 4 5 6 7 8];          % for testing 

SBA.event.condition            = [1:numel(SBA.event.anglesToRotate_e)];         % event conditions
SBA.event.duration              = p.stim.event.duration(2);                     % duration of texture segementation event 
SBA.event.Shapes                = 1:5;                                          % Shape indexes
%SBA.event.ShapeNames            = {'triangle';'square';'diamond';'cross';'L'};
SBA.event.ShapeNames            = {'mask';'dirttruck';'hat';'sumo';'cactus'};   % arbitrary shape names

SBA.event.numShapes            = max(SBA.event.Shapes);                         % one target shape, 4 distractor shapes; 640/2/4/8 = 10 Events per distractor shape
SBA.event.targetShape          = targetshape;                                   % target shape index; 

% fixation cross % NEEDS WORK
p.crs.color                 = [0.5 0.5 0.5 1];          % color of fixation cross
p.crs.size                  = 12;                       % size of fixation cross
p.crs.width                 = 2;                        % width of fixation cross
p.crs.cutout                = 1;                        % are RDK dots allowed to overlapp if fixation cross? 0 = no, 1=yes

% logfiles
if IsLinux; p.log.path              = '/home/pc/matlab/user/sebastian/DUC_Tiltanic/Tiltanic/logfiles/'; 
elseif IsOSX; p.log.path              = '/Users/sebastianwehle/Documents/MATLAB/Data_Tiltanic/logfiles/'; end
p.log.exp_name          = 'SSVEP_Tiltanic';
p.log.add               = '_a';

%% check for logfile being present
filecheck=dir(sprintf('%sVP%02.0f_timing*',p.log.path,p.sub));
if ~isempty(filecheck)
    reply = input(sprintf('\nVP%02.0f existiert bereits. Datei überschreiben? [j/n]... ',p.sub),'s');
    if strcmp(reply,'j')
        p.filename = sprintf('VP%02.0f_timing',p.sub);
    else
        [temp name_ind]=max(cellfun(@(x) numel(x), {filecheck.name}));
        p.filename = sprintf('%s%s',filecheck(name_ind).name(1:end-4),p.log.add);
    end
else
    p.filename = sprintf('VP%02.0f_timing',p.sub);
end

%% Screen initialization
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,0);

% define center of screen
ps.center = [ps.xCenter ps.yCenter];

% make fixation cross textures
[p.FixTex] = MakeFixationTextures(ps, p);
p.crs.rects = CenterRectOnPoint([0 0 p.crs.size p.crs.size], ps.xCenter, ps.yCenter);
Screen('DrawTexture', ps.window, p.FixTex, [], p.crs.rects);
Screen('Flip', ps.window);
%% keyboard and ports setup
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('b') KbName('j') KbName('f') KbName('n')]; % b is unused 
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.TARGET, key.DISTRACTOR, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5),Buttons(6), Buttons(7));
%% shape presentation
if flag_training
    ShapePres.SBA = SBA;
    ShapePres.SBA.event.anglesToRotate = SBA.event.anglesToRotate_t;
    fprintf('####SHAPE PRESENTATION\nFormen werden auf dem Display gezeigt...')
    flag_rep = 1; rep_counter = 0;
    while  flag_rep ==1 
        rep_counter = rep_counter+1;
        ShapePres.presentationOrder = [ShapePres.SBA.event.Shapes(targetshape) ShapePres.SBA.event.Shapes(~ismember(ShapePres.SBA.event.Shapes, targetshape))];
        for iShape = ShapePres.presentationOrder
            ShapePres.SBAin.trial = struct( ...
                    'frames',1,...
                    'cue',1);
            ShapePres.SBAin.trial.event = struct( ...
                    'duration', 1, ...
                    'onset',1,...
                    'contrast',5, ...
                    'direction',randsample([-1 1],1,1),...
                    'shape',iShape); 
            
            ShapePres.Text.Target = '\nDies ist die gesuchte Zielform.\nWeiter mit Leertaste.';
            ShapePres.Text.Distractor = '\nDies ist eine der anderen Formen, die vorkommen können.';
            ShapePres.Text.Rep =  'Wollen Sie die Formen ein weiteres Mal sehen (j/n)?';
        
            if iShape == ShapePres.SBA.event.targetShape
                text2present = ShapePres.Text.Target;
            else 
                text2present = ShapePres.Text.Distractor;
            end 
            % Call the texture generation function
            [barTex, dstRects, angles] = generateBarTextures(ps, ShapePres.SBA, ShapePres.SBAin);
            
            % draw text and background stimuli
            Screen('TextSize', ps.window, 36);
            DrawFormattedText(ps.window, text2present, 'center', [], p.crs.color); % sy = [] at top
            
            Screen('DrawTextures', ps.window, barTex{1}, [], dstRects, angles(:, 1));
            Screen('Flip', ps.window);
        
            inp.prompt_check2 = 0;
            WaitSecs(.2);
            while inp.prompt_check2 == 0             % loop to check for correct input
                [key.keyisdown,key.secs,key.keycode] = KbCheck;
                if key.keycode(key.SPACE)==1
                    inp.prompt_check2 = 1;
                end
                Screen('TextSize', ps.window, 36);
                DrawFormattedText(ps.window, text2present, 'center', [], p.crs.color); % sy = [] at top
        
                Screen('DrawTextures', ps.window, barTex{1}, [], dstRects, angles(:, 1));
                Screen('Flip', ps.window);
            end
        end 
    
        inp.prompt_check1 = 0;
        while inp.prompt_check1 == 0             % loop to check for correct input
            Screen('TextSize', ps.window, 36);
            DrawFormattedText(ps.window, ShapePres.Text.Rep, 'center', 'center', p.crs.color);
            Screen('Flip', ps.window);
            % check for repeated display
            [key.keyisdown,key.secs,key.keycode] = KbCheck;
            if key.keycode(key.TARGET)==1  % yes
               fprintf('\nPressed j – repeating shapes');
               flag_rep = 1; inp.prompt_check1 = 1; 
               fprintf('\nFormen werden wiederholt auf dem Display gezeigt...')
            elseif key.keycode(key.NO)==1
                fprintf('\nPressed n – skipping shapes');
                flag_rep = 0; inp.prompt_check1 = 1;
                fprintf('\nWeiter zu Training...')
            end
            Screen('Flip', ps.window);
        end
        if rep_counter >3
            break;
        end 
    end 
end 
WaitSecs(0.2); % Necessary, otherwise KbCheck assumes previous button press as input

%% inintialize variables
timing = []; button_presses = []; resp = []; randmat = [];
%% training 
if flag_training
    flag_training
    fprintf(1,'\nTraining starten? (j/n)')
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if key.keycode(key.TARGET)==1 % yes
%             fprintf('\n j pressed – starting training1');
            flag_trainend = 0; inp.prompt_check = 1;
        elseif key.keycode(key.NO)==1
%             fprintf('\n n pressed – skipping training1');
            flag_trainend = 1; inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    
    if ~exist('i_tr'); i_tr = 1; end
    while flag_trainend == 0 % do training until ended
        rand('state',p.sub*randi(i_tr*p.sub)) % determine randstate
        [randmat.training, SBA] = rand_Tiltanic(p, RDK, SBA, 1);    % randomization, SBA passed through for training contrasts
        % start experiment
        [timing.training{i_tr},button_presses.training{i_tr},resp.training{i_tr}] = ...
            pres_Tiltanic(p, ps, key, RDK, SBA, randmat.training, i_tr, 1);
        % save logfiles
        save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK', 'SBA')       
        pres_feedback_Tiltanic(resp.training{i_tr},p,ps, key )    

        %loop for training to be repeated
        fprintf(1,'\nTraining wiederholen? (j/n)')
        inp.prompt_check = 0;
        while inp.prompt_check == 0             % loop to check for correct input
            [key.keyisdown,key.secs,key.keycode] = KbCheck;
            if key.keycode(key.TARGET)==1
%                 fprintf('\n j pressed – starting training2');
               i_tr = i_tr + 1; flag_trainend = 0; inp.prompt_check = 1;
            elseif key.keycode(key.NO)==1
%                 fprintf('\n n pressed – skipping training2');
                flag_trainend = 1; inp.prompt_check = 1;
            end
            Screen('Flip', ps.window, 0);
        end
    end
end 
%% experiment: present each block
% randomization
fprintf('\n - - - - - - - - starting experiment  - - - - - - - - \n')
rand('state',p.sub);                                    % determine randstate
[randmat.experiment, SBA] = rand_Tiltanic(p, RDK, SBA, 0);    % randomization, SBA passed through for exp contrasts
for i_tr = 1:randmat.experiment.totaltrials % i_bl = "Trial"
    % start experiment
    [timing.experiment{i_tr},button_presses.experiment{i_tr},resp.experiment{i_tr}] = ...
        pres_Tiltanic(p, ps, key, RDK, SBA, randmat.experiment, i_tr,0);
    % save logfiles
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK', 'SBA')       
    pres_feedback_Tiltanic(resp.experiment{i_tr},p,ps, key )    
end

fprintf(1,'\n\nENDE\n')
ListenChar(0);
sca;

end 