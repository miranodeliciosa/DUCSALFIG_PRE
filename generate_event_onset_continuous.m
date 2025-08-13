function [Event_onsets, conditions] = generate_event_onset_continuous(p, conmat, RDK, SBA, stimFreq, flag_plot)
% generate_event_onset_continuous - Generates temporally jittered event onsets 
% for background and randomly distributed events for foreground stimuli 
% across multiple trials in an SSVEP-style experiment with continous flicker 
% stimulation and trials with a long duration of approx 2 min.
%
% Syntax:
%   [Event_onsets, conditions] = generate_event_onset_continuous( ...
%       nEventsPerCondition, nCond, nContrasts, nTrials, ...
%       frameRate, stimFreq, flag_plot)
%
% Inputs:
%   nEventsPerCondition - Number of background events per contrast condition 
%                         (total across all trials)
%   nCond              - Number of task conditions (e.g., target/non-target)
%   nContrasts         - Number of contrast levels
%   nTrials            - Number of trials to generate
%   frameRate          - Screen refresh rate (in Hz)
%   stimFreq           - [foregroundFreq, backgroundFreq]; flicker frequencies 
%                        for foreground (irrelevant) and background (task-relevant)
%                        If only one value is given, only background events are generated.
%   flag_plot          - Boolean flag; if true, plots and SOA tables will be shown
%
% Outputs:
%   Event_onsets       - Cell array of length nTrials, each containing a struct with:
%       .foreground    - Frame indices of foreground events
%       .background    - Frame indices of background events
%
%   conditions         - Cell array of length nTrials, each containing a struct with:
%       .condition_idx - Condition index (e.g., target = 1, non-target = 2)
%       .contrast_idx  - Assigned contrast level index per background event
%
% Description:
%   This function generates event timings in frame units for both task-relevant 
%   background stimuli and task-irrelevant foreground stimuli in an experiment with 
%   periodic flicker. Background events are spaced with jittered SOAs between 800 ms 
%   and 2000 ms (uniform), while foreground events are randomly distributed across the
%   trial duration while maintaining a minimum SOA.
%
%   Trial durations are extended automatically if 80 background events cannot be 
%   placed with the required SOA constraints. Foreground events are never placed 
%   at the same frames as background events.
%
% Example usage:
%   nEventsPerCondition = 40;
%   nCond = 2;
%   nContrasts = 8;
%   nTrials = 8;
%   frameRate = 60;
%   stimFreq = [7.5, 12];   % Foreground and background flicker frequencies in Hz
%   flag_plot = true;
%
%   [Event_onsets, conditions] = generate_event_onset_continuous( ...
%       nEventsPerCondition, nCond, nContrasts, nTrials, ...
%       frameRate, stimFreq, flag_plot);
%
%   % Access background onsets for trial 3:
%   bg_frames = Event_onsets{3}.background;
%
%   % Access assigned conditions:
%   contrasts = conditions{3}.contrast_idx;
%   shapes = conditions{3}.condition_idx;
%
% @S. Wehle 
%  Leipzig, August 2025

% assign variables from structures
nEventsPerCondition = conmat.nEventsPerCondition;
nCond = numel(p.stim.condition);
nContrasts = numel(p.stim.contrast);
nTrials = conmat.totaltrials;
frameRate = p.scr_refrate;


% Set up constants
minSOA_ms = p.stim.event.SOA_ms(1)*1000;
maxSOA_ms = p.stim.event.SOA_ms(2)*1000;
avgSOA_ms = mean([minSOA_ms maxSOA_ms]);
silence_ms = p.stim.event.min_offset; % 1 second at beginning and end

% Total number of background and foreground events
nEvents_bg_total = conmat.totalSBAevents;
nEvents_fg_total = conmat.totalRDKevents;

nEventsPerTrial_bg = nEvents_bg_total / nTrials;


if mod(nEventsPerTrial_bg,1) ~= 0
    error('Number of background events per trial must be an integer.');
end

% Convert SOA and silence times to frames
minSOA_frames = round(minSOA_ms / 1000 * frameRate);
maxSOA_frames = round(maxSOA_ms / 1000 * frameRate);
silence_frames = round(silence_ms * frameRate);

% Determine flicker steps and possible frames
hasForeground = length(stimFreq) == 2;

if hasForeground
    flickerStep_fg = frameRate / stimFreq(1);
end
flickerStep_bg = frameRate / stimFreq(end);

% Initialize output
Event_onsets = cell(1, nTrials);
conditions = cell(1, nTrials);

% Generate all background condition assignments
all_conditions = repmat(1:nCond, 1, nEventsPerCondition * nContrasts); % see run script to interpret 1 and 2
all_shapes = nan(size(all_conditions));
all_shapes(all_conditions == 1) = SBA.event.targetShape;
nonTargetShape = repmat([SBA.event.Shapes(SBA.event.Shapes ~= SBA.event.targetShape)],1,nEventsPerCondition * nContrasts / numel([SBA.event.Shapes(SBA.event.Shapes ~= SBA.event.targetShape)]));
all_shapes(all_conditions ~= 1) = nonTargetShape;
all_contrasts = repelem(1:nContrasts, nEventsPerCondition * nCond);
all_oridirect = repmat([-1 1], 1, nEventsPerCondition*nContrasts);

idx = randperm(length(all_conditions));
all_conditions = all_conditions(idx);
all_shapes = all_shapes(idx);
all_contrasts = all_contrasts(idx);
all_oridirect = all_oridirect(idx);

% Generate all foreground events (only need number, no condition)
all_fg_event_counts = distribute_events_randomly(nEvents_fg_total, nTrials);
all_directions = Shuffle(repmat([1:4],1, conmat.totalRDKevents/4));
all_eventRDKtype = ones(1,conmat.totalRDKevents);
all_eventRDK = ones(1,conmat.totalRDKevents);

% Create figure for plotting if needed
if flag_plot
    figure; 
    set(gcf, 'Position', [100 100 800 1000])
end

% Main loop over trials
bg_cond_counter = 1;
fg_event_counter = 1;
for trial = 1:nTrials
    trial_ok = false;
    trial_duration_frames = frameRate * p.stim.trialTime; % Initial guess

    while ~trial_ok
        % Compute possible event frames (use same duration for fg and bg)
        possible_bg = round(silence_frames(1)):round(flickerStep_bg):round(trial_duration_frames - silence_frames(2));
        if hasForeground
            possible_fg = round(silence_frames(1)):round(flickerStep_fg):round(trial_duration_frames - silence_frames(2));
        else
            possible_fg = [];
        end

        % Try to place background events with jittered SOAs
        [bg_events, success] = place_jittered_events(possible_bg, nEventsPerTrial_bg, minSOA_frames, maxSOA_frames);

        if success
            % Random placement of foreground events over the full range
            n_fg_trial = all_fg_event_counts(trial);
            fg_events = place_random_events(possible_fg, n_fg_trial, minSOA_frames, bg_events);
            success_fg = ~isempty(fg_events);

            if success_fg
                trial_ok = true;
            else
                trial_duration_frames = trial_duration_frames + round(0.1 * frameRate); % Add 100ms
            end
        else
            trial_duration_frames = trial_duration_frames + round(0.1 * frameRate); % Add 100ms
        end
    end
    
    % save trial duration 
    Event_onsets{trial}.trial_frames = trial_duration_frames;
    Event_onsets{trial}.trial_times = trial_duration_frames/frameRate;


    % Save events
    Event_onsets{trial}.background = bg_events;
    Event_onsets{trial}.background_times = Event_onsets{trial}.background./p.scr_refrate;
    Event_onsets{trial}.foreground = sort(fg_events);
    Event_onsets{trial}.foreground_times = Event_onsets{trial}.foreground/p.scr_refrate;

    % Save conditions for background events
    c_idx = bg_cond_counter:(bg_cond_counter + nEventsPerTrial_bg - 1);
    conditions{trial}.condition_idx = all_conditions(c_idx);
    conditions{trial}.contrast_idx = all_contrasts(c_idx);
    conditions{trial}.shapes_idx = all_shapes(c_idx);
    conditions{trial}.contrOri_idx = all_oridirect(c_idx);
    bg_cond_counter = bg_cond_counter + nEventsPerTrial_bg;

    % Save coherent motion direction for foreground events
    d_idx = fg_event_counter:(fg_event_counter + all_fg_event_counts(trial) -1);
    conditions{trial}.eventRDKtype = all_eventRDKtype(d_idx);
    conditions{trial}.eventRDK = all_eventRDK(d_idx);
    conditions{trial}.direction = all_directions(d_idx);
    fg_event_counter = fg_event_counter + all_fg_event_counts(trial);

    % Plot if requested
    if flag_plot
        subplot(nTrials, 1, trial);
        hold on;
        stem(bg_events / frameRate, ones(size(bg_events)), 'r', 'DisplayName', 'Background');
        stem(fg_events / frameRate, ones(size(fg_events)) * 0.5, 'b', 'DisplayName', 'Foreground');
        legend;
        title(sprintf('Trial %d', trial));
        xlabel('Time (s)');
        ylabel('Event');
        hold off;

        % Print SOA table
        soa = diff(sort(bg_events)) / frameRate * 1000;
        soaTable = table((1:length(soa))', soa', 'VariableNames', {'EventIndex', 'SOA_ms'});
        fprintf('Trial %d Background SOAs (ms):\n', trial);
        disp(soaTable);
    end
end
end

%% HELPER FUNCTIONS %%%%%%

% --- Helper function to place jittered events ---
function [events, success] = place_jittered_events(possibleFrames, nEvents, minSOA, maxSOA, excludeFrames)
    if nargin < 5
        excludeFrames = [];
    end

    events = [];
    remainingFrames = setdiff(possibleFrames, excludeFrames);
    remainingFrames = sort(remainingFrames);
    
    if length(remainingFrames) < nEvents
        success = false;
        return;
    end

    currentFrame = remainingFrames(1);
    events = currentFrame;
    for i = 2:nEvents
        minNext = currentFrame + minSOA;
        maxNext = currentFrame + maxSOA;
        candidates = remainingFrames(remainingFrames >= minNext & remainingFrames <= maxNext);
        if isempty(candidates)
            success = false;
            return;
        end
        currentFrame = candidates(randi(length(candidates)));
        events(end+1) = currentFrame;
    end
    success = true;
end

% --- Helper function to place random foreground events ---
function events = place_random_events(possibleFrames, nEvents, minSOA, excludeFrames)
    events = [];
    remainingFrames = setdiff(possibleFrames, excludeFrames);
    remainingFrames = sort(remainingFrames);

    if length(remainingFrames) < nEvents
        return;
    end

    % Randomly permute and filter with minimum SOA
    permFrames = remainingFrames(randperm(length(remainingFrames)));
    events = permFrames(1);

    for i = 2:length(permFrames)
        if all(abs(permFrames(i) - events) >= minSOA)
            events(end+1) = permFrames(i);
            if length(events) == nEvents
                break;
            end
        end
    end

    if length(events) < nEvents
        events = [];
    end
end

% --- Helper to distribute foreground events randomly across trials ---
function trial_counts = distribute_events_randomly(total_events, nTrials)
    trial_counts = zeros(1, nTrials);
    for i = 1:total_events
        idx = randi(nTrials);
        trial_counts(idx) = trial_counts(idx) + 1;
    end
end
