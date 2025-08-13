function [ conmat ] = rand_DuCSalFig(p, RDK, SBA, flag_training)

% set trial number and check whether disribution works
if flag_training~=0
    conmat.totaltrials = numel(p.stim.condition)*numel(p.stim.eventnum)*p.stim.con_repeats_t;
    conmat.totalblocks = 1;
    if mod(sum(p.stim.eventnum)*p.stim.con_repeats_t,numel([1 p.stim.RDKdistr]))~=0
        error('rando:targdistrdistribute', ['Can not distribute target and distractor types equally across events. The event number per condition is %1.0f with %1.0f repetition(s).\n',...
        '%1.0f types need to be distributed. Recommended to increase con_repeats to e.g. %1.0f.'],... 
        sum(p.stim.eventnum), p.stim.con_repeats_t, numel([1 p.stim.RDKdistr]), lcm(sum(p.stim.eventnum),numel([1 p.stim.RDKdistr]))/sum(p.stim.eventnum));
    end
    conmat.nEventsPerCondition=p.stim.con_repeats_t;
else
    conmat.totalSBAevents = numel(p.stim.contrast)*numel(p.stim.condition)*p.stim.con_repeats_e;
    conmat.totalRDKevents = p.stim.event.ratio*conmat.totalSBAevents;
    conmat.totaltrials = conmat.totalSBAevents/p.stim.eventsPerTrial;
    conmat.nEventsPerCondition=p.stim.con_repeats_e;
%     if mod(sum(p.stim.eventnum)*p.stim.con_repeats_e,numel([1 p.stim.RDKdistr]))~=0
%         error('rando:targdistrdistribute', ['Can not distribute target and distractor types equally across events. The event number per condition is %1.0f with %1.0f repetition(s).\n',...
%         '%1.0f types need to be distributed. Recommended to increase con_repeats to e.g. %1.0f.'],... 
%         sum(p.stim.eventnum), p.stim.con_repeats_e, numel([1 p.stim.RDKdistr]), lcm(sum(p.stim.eventnum),numel([1 p.stim.RDKdistr]))/sum(p.stim.eventnum));
%     end   

end 

% matrix with onset times of on frames for RDKs
t.onframesonset = nan(numel(RDK.RDK),p.scr_refrate*p.stim.trialTime);
t.onframesonset_times = t.onframesonset; % onset times in s
for i_rdk = 1:numel(RDK.RDK)    
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2)); %MD: bei welchen frames geht ein RDK "an"; für alle RDKs/freqs; edit: bei Frequenzen bei denen es nicht aufgeht: zum nächsten Frame aufgerundet.
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end
% brownian motion movement of RDK dots
t.movonset_frames=nan(1,p.scr_refrate*p.stim.trialTime); % brownian motion at all frames continously 
t.movonset_times=nan(1,p.scr_refrate*p.stim.trialTime);
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;

%% generate randomized trial parameters for background and foreground stimuli 
flag_plot = 0;

% event_onsets{trials}.background and foreground fields 
[event_onsets, conmat.mats] = generate_event_onset_continuous(p, conmat, RDK, SBA, [RDK.RDK.freq SBA.freq], flag_plot);

%% write all information into trial structure

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
    conmat.trials(i_tr).SBA_eventshapeNames = SBA.event.ShapeNames{conmat.mats{i_tr}.shapes_idx};

end

end