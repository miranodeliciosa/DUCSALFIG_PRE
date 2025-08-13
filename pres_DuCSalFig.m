function [timing,key,resp] = pres_DuCSalFig(p, ps, key, RDK, SBA, conmat, i_tr, flag_training)
% presents experiment SSVEP_DuCSalFig
%   p               = parameters
%   ps              = screen parameters
%   RDK             = RDK parameters
%   SBA             = SBA parameters (Static Bar Array)
%   blocknum        = number of block
%   flag_training   = flag for trainig (1) or experiment (0)

%% adaptations for training
if flag_training == 1
    blocknum_present = blocknum;
    blocknum = 1;
end

%% initialize some RDK settings
% define input for RDK init function
RDKin.scr = ps; 
RDKin.scr.refrate = p.scr_refrate;
RDKin.Propixx = p.scr_imgmultipl;
RDKin.RDK = RDK;
RDKin.crs = p.crs;

%% loop for each trial
% send start trigger
if flag_training~=1
%     PRPX_sendRecTrigger('start') % TODO: how can i send trigger in 119?
elseif flag_training == 1 && p.stim.train_trials ~= 0
    trialindex = [1:p.stim.train_trials]'; %quick fix to make training shorter; 
end

% Wait for release of all keys on keyboard, then sync us to retrace:
KbWait(ps.RespDev,1); % wait for releasing keys on indicated response device

% create keayboard queue
KbQueueCreate(ps.RespDev)


% Build the German feedback text
text2present = [ ...
    '\nR E A D Y ?' ...
    '\n\nDrücken Sie die Leertaste, um fortzufahren.'
];

% Draw
Screen('TextSize', ps.window, 24);
DrawFormattedText(ps.window, text2present, 'center', 'center', p.crs.color);

% Show it
Screen('Flip', ps.window);

% [key.pressed, key.firstPress]=KbQueueCheck;
key.rkey=key.SPACE;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    DrawFormattedText(ps.window, text2present, 'center', 'center', p.crs.color);
    Screen('Flip', ps.window, 0);                   % flip screen
end

Screen('Flip', ps.window, 0); 
%experimenter output
if flag_training~=1
    fprintf(1,'\nexperiment - Praesentation - Trial %1.0f', i_tr)
else
    fprintf(1,'\ntraining - Praesentation - Trial %1.0f', i_tr)
end

ttt=WaitSecs(1);

% loop across trials

   % fprintf('%1.0f',mod(i_tr,10))
    inittime=GetSecs;
    %% initialize trial structure, RDK, SBA, cross, logs
    %RDK
    RDKin.trial = struct('duration',conmat.trials(i_tr).trial_time,...
        'frames',conmat.trials(i_tr).trial_frames,...
        'cue',1); % start of trial?
    RDKin.trial.event = struct('onset',conmat.trials(i_tr).event_onset_frames,...
        'direction',conmat.trials(i_tr).eventdirection,...
        'RDK',conmat.trials(i_tr).eventRDK);

    [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init_FShiftGlob(RDKin.scr,RDKin.Propixx,RDKin.RDK,RDKin.trial,RDKin.crs);
    
    %SBA
    SBAin.trial = struct( ...
        'frames',RDKin.trial.frames,...
        'cue',1);

    SBAin.trial.event = struct( ...
        'duration', SBA.event.duration, ...
        'onset',conmat.trials(i_tr).SBA_event_frames,...%!!!
        'contrast',conmat.trials(i_tr).SBA_contrast, ...
        'direction', conmat.trials(i_tr).SBA_contrast_direction,...
        'shape',conmat.trials(i_tr).SBA_eventshapes); % EXCLUDE NANs

    [barTex, dstRects, angles] = generateBarTextures(ps, SBA, SBAin);

    % crossmat omitted 
    % preallocate timing
    timing(i_tr) = struct('VBLTimestamp',NaN(1,frames.flips),'StimulusOnsetTime',NaN(1,frames.flips),...
        'FlipTimestamp',NaN(1,frames.flips),'Missed',NaN(1,frames.flips));

%% set up responses
    %setup key presses
    key.presses{1}=nan(size(colmat,3),sum(key.keymap));
    key.presses_t{1}=nan(size(colmat,3),sum(key.keymap));
    
    % TODO set up responses for SBA events and documentation
    resp(1).trialnumber                  = i_tr;
    resp(1).trial_time                   = conmat.trials(i_tr).trial_time;
    resp(1).trial_frames                 = conmat.trials(i_tr).trial_frames;
    resp(1).eventtype                    = conmat.trials(i_tr).eventtype;
    resp(1).eventRDK                     = conmat.trials(i_tr).eventRDK;
    resp(1).eventdirection               = conmat.trials(i_tr).eventdirection;
    resp(1).RDK_event_onset_frames       = conmat.trials(i_tr).event_onset_frames;
    resp(1).RDK_event_onset_times        = conmat.trials(i_tr).event_onset_times;
    resp(1).SBA_event_frames             = conmat.trials(i_tr).SBA_event_frames;
    resp(1).SBA_event_times              = conmat.trials(i_tr).SBA_event_times;
    resp(1).SBA_contrast                 = conmat.trials(i_tr).SBA_contrast;
    resp(1).SBA_eventtype                = conmat.trials(i_tr).SBA_eventtype;
    resp(1).SBA_eventshapes              = conmat.trials(i_tr).SBA_eventshapes;
    resp(1).SBA_eventshapenames          = {conmat.trials(i_tr).SBA_eventshapeNames};

    %% set up trigger vector
    % TODO...

     % draw fixation cross
     Screen('DrawTexture', ps.window, p.FixTex, [], p.crs.rects);
     %% keyboard
    % start listening to keyboard
    KbQueueStart(ps.RespDev);
    KbQueueFlush(ps.RespDev);
    
    % flip to get everything synced
    Screen('Flip', ps.window, 0);

       %% loop across flips   
    for i_fl = 1:frames.flips
        %% Drawing

        % SBA - background 
        Screen('DrawTextures', ps.window, barTex{i_fl}, [], dstRects, angles(:, i_fl));
 
        % RDK - foreground
        Screen('DrawDots', ps.window, dotmat(:,:,i_fl), dotsize(:,i_fl), colmat(:,:,i_fl), ps.center, 0, 0);

        

        % fixation cross
        Screen('DrawTextures', ps.window, p.FixTex, [], p.crs.rects);
        
        %% start trigger schedule and start listening to response device
        if i_fl == 1 % send the trigger with the start of the 1st flip
%             Datapixx('RegWrVideoSync');
        end
        
        % Flip
        [timing(i_tr).VBLTimestamp(i_fl), timing(i_tr).StimulusOnsetTime(i_fl), timing(i_tr).FlipTimestamp(i_fl), timing(i_tr).Missed(i_fl)] = Screen('Flip', ps.window, 0);
        
        % save image
        %imageArray=Screen('GetImage', ps.window);
        %imwrite(imageArray,'ScreenTest_cue_FShiftGlob.png')
        
        % send trigger/save timing/ reset timing
        if i_fl == 1
            % start trigger
            starttime=GetSecs;
            KbEventFlush(ps.RespDev); % flush keyboard
        end
        
        %% check for button presses
        [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
        key.presses{1}(i_fl,:)=key.firstPress(key.keymap)>1;
        key.presses_t{1}(i_fl,:)=(key.firstPress(key.keymap)-starttime).*key.presses{1}(i_fl,:);
%         if any(key.firstPress(key.keymap)>1)
%             lptwrite(1,find(key.firstPress(key.keymap),1,'first'),500);
%         end
    % log flips
    end 
     %% ITI
    % draw fixation cross again
    Screen('DrawTextures', ps.window, p.FixTex, [], p.crs.rects);

    % flip
    Screen('Flip', ps.window, 0);
    
    % get time
    crttime = GetSecs;
    
    % add waiting period to check for late button presses -  already
    % accounted for in rando
    %ttt=WaitSecs(p.targ_respwin(2)/1000-p.stim.event.min_offset(2)-p.stim.event.length);
    
    % check for button presses
    [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
    key.presses{1}(i_fl+1,:)=key.firstPress(key.keymap)>1;
    key.presses_t{1}(i_fl+1,:)=(key.firstPress(key.keymap)-starttime).*key.presses{1}(i_fl+1,:);
    
    % do behavioral calculation
    % for CODING
% % %     iTr = 1;
% % %     key.keymap_ind = button_presses.experiment{1}.keymap_ind;
% % %     key.TARGET = button_presses.experiment{1}.TARGET;
% % %     key.DISTRACTOR = button_presses.experiment{1}.DISTRACTOR;
% % %     key.presses = button_presses.experiment{iTr}.presses;
% % %     key.presses_t = button_presses.experiment{iTr}.presses_t;
% % % 
% % %     resp.SBA_event_times = resp(1).experiment{iTr}.SBA_event_times;
% % %     resp.SBA_eventtype = resp(1).experiment{iTr}.SBA_eventtype;
% % %    

    key.presses{1}(1,:)=[]; % erase first row?
    key.presses_t{1}(1,:)=[]; 
    
    % get frame and timing of button press onsets
    resp(1).button_presses_fr=nan(max(sum(key.presses{1})),size(key.presses{1},2));
    resp(1).button_presses_t=resp(1).button_presses_fr;
    for i_bt = 1:size(key.presses{1},2)
        try
            resp(1).button_presses_fr(1:sum(key.presses{1}(:,i_bt)),i_bt)=...
                find(key.presses{1}(:,i_bt));
            resp(1).button_presses_t(1:sum(key.presses{1}(:,i_bt)),i_bt)=...
                key.presses_t{1}(find(key.presses{1}(:,i_bt)),i_bt)*1000; % in ms
        catch
            resp(1).button_presses_fr(:,i_bt)=nan;
        end
    end    

     % check for hits {'hit','miss','CR','FA_proper','FA'} %ToDo: introduce more FA types??
    resp(1).button_presses_type = {}; %{'hit','miss','CR','FA_proper','FA'}
    resp(1).button_presses_class = {}; %{'correct';'incorrect'}
    resp(1).button_presses_RT = []; % reaction time in ms (after event or closest to other event)
    resp(1).event_response_class = {}; %{'correct';'incorrect'}
    resp(1).event_response_type = {}; %{'hit','miss','CR','FA_proper'}
    resp(1).event_response_RT = []; %reaction time or nan
    resp(1).button_pressed = [];
    resp(1).button_pressed_event = [];

    % all relevant presses 
    t.presses = resp(1).button_presses_t(:,key.keymap_ind==key.TARGET | key.keymap_ind==key.DISTRACTOR);
    % first define response windows
    t.respwin = (resp(1).SBA_event_times'+(p.targ_respwin/1000))*1000; % reponse times in ms!
    t.postcue_presses = 0; %index presses

    % ---- Prep: event onsets (ms), types, windows ----
    temp.e_on_ms   = resp(1).SBA_event_times(:) * 1000;        % [nEvents x 1], ms
    temp.e_type    = resp(1).SBA_eventtype(:);                 % [nEvents x 1], 1=target,2=distractor
    temp.numEvents = numel(resp(1).SBA_eventtype(:));
    
    % ---- Prep: presses flattened but we write results back in the same shape as t.presses ----
    [temp.maxP, temp.nButtons] = size(t.presses);                    % nButtons should be 2 (1='k', 2='n')
    
    % Precompute button codes per column for writing back in-place
    temp.btnCodesMat = nan(temp.maxP, temp.nButtons);
    if temp.nButtons >= 1, temp.btnCodesMat(:,1) = 1; end           % 1 = target button ('k')
    if temp.nButtons >= 2, temp.btnCodesMat(:,2) = 2; end           % 2 = distractor button ('n')
    
    % Flatten valid presses (chronological tie-breaker across buttons)
    temp.all_idx_lin  = find(~isnan(t.presses));
    temp.all_times    = t.presses(temp.all_idx_lin);
    [~, temp.ord]     = sort(temp.all_times, 'ascend');
    temp.lin_ordered  = temp.all_idx_lin(temp.ord);
    temp.times_ordered= temp.all_times(temp.ord);
    temp.btn_ordered  = temp.btnCodesMat(temp.lin_ordered);
    
    % ---- Outputs (initialize) ----
    resp(1).button_presses_type  = repmat({NaN}, temp.maxP, temp.nButtons); % per-press label in t.presses order
    resp(1).button_presses_class = repmat({NaN}, temp.maxP, temp.nButtons); % per-press "correct"/"incorrect"
    resp(1).button_presses_RT    = nan(temp.maxP, temp.nButtons);           % per-press RT (ms)
    resp(1).button_pressed       = nan(temp.maxP, temp.nButtons);           % per-press button code (1/2)
    
    resp(1).event_response_type  = repmat({NaN}, temp.numEvents , 1);          % per-event label ('hit','miss','CR','FA'), NaN if none
    resp(1).event_response_class = repmat({NaN}, temp.numEvents , 1);          % per-event "correct"/"incorrect", NaN if none
    resp(1).event_response_RT    = nan(temp.numEvents, 1);                    % per-event RT (ms), NaN if none
    resp(1).button_pressed_event = nan(temp.numEvents, 1);                    % per-event button code (1/2), NaN if none
    
    % Track whether an event already got its first press
    temp.event_claimed = false(temp.numEvents,1);
    
    % ---- Main loop over presses in chronological order ----
    for k = 1:numel(temp.lin_ordered)
        lin   = temp.lin_ordered(k);
        temp.t_ms  = temp.times_ordered(k);
        temp.bcode = temp.btn_ordered(k);
    
        % Candidate events whose response window contains this press
        temp.inwin = (temp.t_ms >= t.respwin(:,1)) & (temp.t_ms <= t.respwin(:,2));
    
        % If none, it's a proper false alarm (outside all windows)
        if ~any(temp.inwin)
            resp(1).button_presses_type{lin}  = 'FA_proper';
            resp(1).button_presses_class{lin} = 'incorrect';
            resp(1).button_presses_RT(lin)    = NaN;      % no event, so no RT
            resp(1).button_pressed(lin)       = temp.bcode;
        else % there is a button press
            % If multiple events overlap, choose the one with nearest onset
            temp.candIdx = find(temp.inwin);
            [~, temp.imin] = min(abs(temp.t_ms - temp.e_on_ms(temp.candIdx)));
            ev = temp.candIdx(temp.imin);
        
            % Enforce "first press counts": ignore/link? -> spec says at most one per event;
            % subsequent presses during the same window are treated as FA_proper (spurious).
            if temp.event_claimed(ev)
                resp(1).button_presses_type{lin}  = 'FA_proper';
                resp(1).button_presses_class{lin} = 'incorrect';
                resp(1).button_presses_RT(lin)    = NaN;
                resp(1).button_pressed(lin)       = temp.bcode;            
            end
        
           % ---- Switch–case classification (no labelPress) ----
            % Default (defensive)
            temp.rtype  = 'FA';          % will be overwritten below
            temp.rclass = 'incorrect';
        
            switch temp.e_type(ev)
                case 1   % target event
                    switch temp.bcode
                        case 1
                            temp.rtype  = 'hit';
                            temp.rclass = 'correct';
                        case 2
                            temp.rtype  = 'miss';
                            temp.rclass = 'incorrect';
                        otherwise
                            % unknown button; keep defaults
                    end
                case 2   % distractor event
                    switch temp.bcode
                        case 2
                            temp.rtype  = 'CR';
                            temp.rclass = 'correct';
                        case 1
                            temp.rtype  = 'FA';
                            temp.rclass = 'incorrect';
                        otherwise
                            % unknown button; keep defaults
                    end
                otherwise
                    % unknown event type; keep defaults
            end
            
            % Per-press outputs (aligned to t.presses position)
            resp(1).button_presses_type{lin}  = temp.rtype;
            resp(1).button_presses_class{lin} = temp.rclass;
            resp(1).button_presses_RT(lin)    = temp.t_ms - temp.e_on_ms(ev);
            resp(1).button_pressed(lin)       = temp.bcode;
        
            % Per-event outputs (only first press wins)
            resp(1).event_response_type{ev}   = temp.rtype;
            resp(1).event_response_class{ev}  = temp.rclass;
            resp(1).event_response_RT(ev)     = temp.t_ms - temp.e_on_ms(ev);
            resp(1).button_pressed_event(ev)  = temp.bcode;
        
            event_claimed(ev) = true;
        end 
    end

    % wait
    crttime2 = GetSecs;
    ttt=WaitSecs(.1);
end 