function [ conmat, SBA] = rand_Tiltanic(p, RDK, SBA, flag_training)
% rand_Tiltanic - Generates randomized trial structure for SBA and RDK events
%
%   [conmat, SBA] = rand_Tiltanic(p, RDK, SBA, flag_training)
%
%   This function prepares the full trial matrix for the Tiltanic experiment by assigning
%   randomized timing and conditions for both:
%       • Task-relevant background segmentation events (SBA)
%       • Task-irrelevant foreground RDK motion events
%
%   It adapts the number of trials, contrasts, and event repetitions depending on
%   whether the experiment is run in training or full mode. It also computes flicker
%   timing for both stimulus layers and populates the per-trial data structure with
%   all relevant event information.
%
%   INPUTS:
%       p              - Experiment parameter structure (subject ID, design, screen setup)
%       RDK            - Foreground dot stimulus structure (e.g., motion parameters, frequency)
%       SBA            - Background bar stimulus structure (e.g., size, color, contrast angles)
%       flag_training  - 1 = training mode (fewer events, limited contrast range)
%                        0 = full experiment mode
%
%   OUTPUTS:
%       conmat         - Structure containing trial-wise event assignments and timing:
%                           • Event onsets (frame/time) for background and foreground
%                           • Assigned contrast levels, conditions, shapes, and directions
%                           • Shape names and trial durations
%       SBA            - Updated SBA structure with trial-specific anglesToRotate
%
%   EVENT TIMING:
%       Background and foreground event onsets are computed using
%       `generate_event_onset_continuous`, which ensures:
%           • SOA jitter for background events (1000–4000 ms)
%           • Foreground events are randomly distributed
%           • Temporal separation between layers (no event overlap)
%
%   USAGE EXAMPLE:
%       [conmat, SBA] = rand_Tiltanic(p, RDK, SBA, 1);  % Generate training trials
%
%   DEPENDENCIES:
%       • generateBarTextures.m: Used indirectly during trial preparation for SBA event shapes
%
%   SEE ALSO:
%       run_Tiltanic, pres_DuCSalFig, generateBarTextures
%
%   Author: Sebastian Wehle, Leipzig (2025)


%% Determine Experiment Mode and Set Trial Parameters
% Depending on whether the function is called for training or experiment,
% adjust contrast levels, number of events per condition, and total trial time.

if flag_training~=0
    SBA.event.anglesToRotate       = SBA.event.anglesToRotate_t;
    conmat.totalSBAevents = numel(p.stim.contrast)*numel(p.stim.condition)*p.stim.con_repeats_t;
    conmat.totalRDKevents = p.stim.event.ratio*conmat.totalSBAevents;
    conmat.totaltrials = conmat.totalSBAevents/p.stim.eventsPerTrial_t;
    conmat.nEventsPerCondition=p.stim.con_repeats_t;
    p.stim.trialTime = p.stim.trialTime_t;
else
    SBA.event.anglesToRotate       = SBA.event.anglesToRotate_e;
    conmat.totalSBAevents = numel(p.stim.contrast)*numel(p.stim.condition)*p.stim.con_repeats_e;
    conmat.totalRDKevents = p.stim.event.ratio*conmat.totalSBAevents;
    conmat.totaltrials = conmat.totalSBAevents/p.stim.eventsPerTrial_e;
    conmat.nEventsPerCondition=p.stim.con_repeats_e;
    p.stim.trialTime = p.stim.trialTime_e;
end 
% Update SBA Parameters According to Experiment Mode
SBA = SBA;

%% Generate RDK Flicker Timing (On-frame Onsets)
% Compute the frames at which each RDK dot turns "on", based on flicker frequency.
% Store both frame indices and their corresponding time values.

t.onframesonset = nan(numel(RDK.RDK),p.scr_refrate*p.stim.trialTime);
t.onframesonset_times = t.onframesonset; % onset times in s
for i_rdk = 1:numel(RDK.RDK)    
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2)); %MD: bei welchen frames geht ein RDK "an"; für alle RDKs/freqs; edit: bei Frequenzen bei denen es nicht aufgeht: zum nächsten Frame aufgerundet.
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end
%% Generate RDK Motion Timing (Brownian Movement)
% Determine at which frames the RDK dot positions should update (i.e., move).
% This is independent from flicker and can vary in frequency.

t.movonset_frames=nan(1,p.scr_refrate*p.stim.trialTime); % brownian motion at all frames continously 
t.movonset_times=nan(1,p.scr_refrate*p.stim.trialTime);
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;

%% Generate Randomized Event Onsets for Foreground and Background
% Call generate_event_onset_continuous to produce time-jittered event timings per trial
% for both foreground RDK and background SBA segmentation events.
% Ensure events do not overlap and follow frequency constraints.

[event_onsets, conmat.mats] = generate_event_onset_continuous(p, conmat, RDK, SBA, [RDK.RDK.freq SBA.freq], 0);

%% Build Trial-wise Structure for All Events
% Loop over trials to fill in `conmat.trials(i)` fields:
% timing, shapes, contrasts, directions, event onsets, types, etc.
% Each trial contains both SBA (shape) and RDK (motion) event info. all information into trial structure

for i_tr = 1:conmat.totaltrials
    % trialnumber
    conmat.trials(i_tr).trialnum = i_tr;

    % post-cue times
    conmat.trials(i_tr).trial_time = event_onsets{i_tr}.trial_times;
    
    % post-cue frames
    conmat.trials(i_tr).trial_frames = event_onsets{i_tr}.trial_frames;
    
    % type of events, here only 1
    conmat.trials(i_tr).eventtype = conmat.mats{i_tr}.eventRDKtype;
    
    % which RDK shows event? here only 1
    conmat.trials(i_tr).eventRDK = conmat.mats{i_tr}.eventRDK;
    
    % eventdirection ((according to RDK.event.direction) [1 2 3 4])
    conmat.trials(i_tr).eventdirection = conmat.mats{i_tr}.direction;
    
    % event onset frames foreground RDK
    conmat.trials(i_tr).event_onset_frames = event_onsets{i_tr}.foreground;
    
    % event onset times foreground RDK
    conmat.trials(i_tr).event_onset_times = event_onsets{i_tr}.foreground_times;
    
    % SBA event frames backgorund 
    conmat.trials(i_tr).SBA_event_frames = event_onsets{i_tr}.background;
    
    % SBA event times background 
    conmat.trials(i_tr).SBA_event_times = event_onsets{i_tr}.background_times;
 
    % SBA event contrast 
    conmat.trials(i_tr).SBA_contrast = conmat.mats{i_tr}.contrast_idx;

    % SBA event orientation change direction 
    conmat.trials(i_tr).SBA_contrast_direction = conmat.mats{i_tr}.contrOri_idx;
    
    % SBA event type (target, non-target; [1 2])
    conmat.trials(i_tr).SBA_eventtype = conmat.mats{i_tr}.condition_idx;

    % SBA event shape (target shape [1], non-target [2 3 4 5] )
    conmat.trials(i_tr).SBA_eventshapes = conmat.mats{i_tr}.shapes_idx;

    % SBA event shape name ('triangle';'square';'diamond';'bar hor';'bar vert' )
    conmat.trials(i_tr).SBA_eventshapeNames = {SBA.event.ShapeNames{conmat.mats{i_tr}.shapes_idx}};

end

end