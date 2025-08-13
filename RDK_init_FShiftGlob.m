function [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init_FShiftGlob(scr,Propixx,RDK,trial,cross)
%RDK_init initializes parameters for RDKs
% This function initializes dot colors and dot positions for every frame.
% The different RDKs are treated as "one big dot cloud". Dots are
% randomized, to ensure that dots of lastly drawn RDK are not systematically
% overlapping the others.
% requires
%       scr.        refrate
%                   shift           shift to 4 quadrants (specific to Propixx)
%       Propixx     4               specific fro Propixx to draw in QUAD4X mode
%       RDK.RDK(x).                 input structure ith RDKs
%                   size            [x y] of ellipse in pixel
%                   centershift     displacement from center in pixel
%                   col             [R G B A; R G B A]; on off colors
%                   mov_freq        frequency for which movement is updated
%                   mov_dir         movement direction [x y; x y;...] in pixel
%                   mov_speed       movement speed in pixel
%                   num             number of dots for each RDK
%                   dot_size        size of dots
%       RDK.event.  type            'globalmotion'; only works for global motion so far
%                   duration        duration in s
%                   coherence       amount of dots showing coherent motion
%                   direction       direction of event movements [x y; x y;...] in pixel
%       trial.      duration        in s
%                   frames
%                   cue             cue frame (time point for syncing)
%                   event.onset     onset frame of event
%                   event.direction direction of event movements [1 2 3 4] in pixel
%                   event.RDK       RDK for event
%       cross.      cutout          no movements close to cross? 1 = no movements
%                   size            size of central fixation cross
   
% ToDo: Put this elsewhere when setting up a trial?
if ~isfield(trial, 'frames')
    frames.pertrial = trial.duration*scr.refrate;
else
    frames.pertrial = trial.frames;
end
frames.flips = frames.pertrial/Propixx; % There are 4 frames per flip in the QUAD4X mode

colmat  = [];
dotmat  = [];
dotsize = [];
rdkidx  = [];

%% Do computations for every RDK

for r = 1:length(RDK.RDK)
    
    onframes = scr.refrate/RDK.RDK(r).freq/2; %number of "on"-frames
    
    % setup color matrix sync on cue
    if ~isfield(trial, 'cue') % check for at certain frame (shifting design) vs cued befor trial design
        fr_cue = 1;
    else
        if trial.cue < 2; fr_cue = 1;
        else
            fr_cue = trial.cue;
        end
    end
    
    % find number of decimals for frequency to create vector that ends with full cycle (for frequencies with decimal point)
    if mod(RDK.RDK(r).freq,1) > 0
        t.temp = strsplit(num2str(RDK.RDK(r).freq),'.');
        t.dec = length(t.temp{2});
    else
        t.dec = 0;
    end
    % create vector of x*2*10^t.dec seconds that contains two parts ending with a full cycle and that it is at least twice as long
    % as the pre-cue or post-cue window
    if any(scr.refrate*10^t.dec < [frames.pertrial-fr_cue+1 fr_cue-1])
        t.framenum = 2 * scr.refrate * 10^t.dec * ...
            max(ceil([frames.pertrial-fr_cue+1 fr_cue-1]/scr.refrate*10^t.dec));
    else
        t.framenum = 2*scr.refrate*10^t.dec;
    end
    % create frame vector for MakeFlicker function
    t.frames = 0:t.framenum-1;
    % create raw weights with interpolated flanks (is longer than it should be
    weight_raw = MakeFlicker(t.frames,scr.refrate/RDK.RDK(r).freq);
    % figure; plot(t.frames./scr.refrate, weight_raw)
    % prune weights according to number "post-cue" frames
%     t.idx = [ (1:fr_cue-1) + numel(weight_raw)/2 - (fr_cue-1)... % pre-cue
%         (numel(weight_raw)/2)+(1:frames.pertrial-fr_cue+1)]; % post-cue
    t.idx = [ (numel(weight_raw)/2)+(1:frames.pertrial)]; % post-cue
    weight = weight_raw(t.idx);
    % figure; plot(t.frames(t.idx)./scr.refrate, weight)
    % figure; plot( weight)
     
%     % compute the colour weights (1 --> 100% color 1; 0 --> 100% color 2)
%      weight_postcue = MakeFlicker(0:frames.pertrial-fr_cue,scr.refrate/RDK.RDK(r).freq); % make sure flicker starts with an on cycle when cue changes
%      weight_precue = MakeFlicker(0:fr_cue-2,scr.refrate/RDK.RDK(r).freq); 
%      weight_precue = flip(abs(weight_precue-1)); % make sure there is a full off cycle before cue change %ToDo: is this flipping ok?
%      weight = [weight_precue weight_postcue];
     
     colmat_r = zeros(size(RDK.RDK(r).col,2), RDK.RDK(r).num, frames.pertrial); %setup color matrix: RGB x dot number x frames
     for w=1:length(weight)
         colmat_r(:,:,w) = repmat((weight(w)*RDK.RDK(r).col(1,:) + (1-weight(w))*RDK.RDK(r).col(2,:))',[1,RDK.RDK(r).num]); %Formula to compute weighted colour for every frame
     end
     
%      temp = colmat_r;
     
%      %ToDo: get rid of this, this is old
%     % post cue
%     colidx_postc = repmat([ones(1,onframes) ones(1,onframes)*2], RDK.RDK(r).num, ceil(numel(fr_cue:frames.pertrial)/(onframes*2))); % write 1 for "on"-frames and 2 for "off"-frames
%     colidx_postc = colidx_postc(:,1:numel(fr_cue:frames.pertrial)); % cut to trial length
%     % pre cue
%     colidx_prec = repmat([ones(1,onframes) ones(1,onframes)*2], RDK.RDK(r).num, ceil(numel(1:fr_cue-1)/(onframes*2))); % write 1 for "on"-frames and 2 for "off"-frames
%     colidx_prec = colidx_prec(:,end-fr_cue+2:end); % cut to trial length
%     % both together
%     colidx = [colidx_prec colidx_postc];
%     % figure; plot(colidx(1,:)); hold on; plot([fr_cue fr_cue], [0.9 2.1])
% 
%    
%     colmat_r = zeros(size(RDK.RDK(r).col,2), RDK.RDK(r).num, frames.pertrial); %setup color matrix: RGB x dot number x frames
%     colmat_r = reshape(RDK.RDK(r).col(colidx,:)',size(colmat_r)); %fill color matrix: write in the RGB triplets for every frame
%     %ToDo: Add here a modulator for interpolated flicker frequencies or sinuoid flickers...
%     %figure;plot(squeeze(temp(1,1,:)));hold on; plot(squeeze(colmat_r(1,1,:)));hold on; plot([fr_cue fr_cue], [0.05 0.16];
%     
    
    % movdot says when to do dot position updates
    if RDK.RDK(r).mov_freq == 0
        movdot = double([abs(diff([2 weight]))==~0]); % Writes a "1" for every first on-frame and every first off-frame
        % figure; plot(colidx(1,:)); hold on; plot([fr_cue fr_cue], [-0.1 2.1]); plot(movdot)
    else
        movdot_update = scr.refrate/RDK.RDK(r).mov_freq;
        movdot = zeros(1,size(weight,2));
        movdot(fr_cue:movdot_update:end)=1; % again synced at cue onset
        movdot(fr_cue:-movdot_update:1)=1;
        % figure; plot(colidx(1,:)); hold on; plot([fr_cue fr_cue], [-0.1 2.1]); plot(movdot)
    end
    
    % initalize RDK positions
    RDK.RDK(r).rad_x = RDK.RDK(r).size(1)/2-RDK.RDK(r).dot_size/2; % Define the maximal x and y positions for a dot
    RDK.RDK(r).rad_y = RDK.RDK(r).size(2)/2-RDK.RDK(r).dot_size/2;
    cross.rad = (cross.size/2+RDK.RDK(r).dot_size)*cross.cutout; % Define the minimal x and y positions for a dot (cutout around fixation cross) %ToDo: is it ok to calculate with inf? cause cross.rad may be 0 and is a divisor later
    dot_pos = [randi([-RDK.RDK(r).rad_x RDK.RDK(r).rad_x],RDK.RDK(r).num,1) randi([-RDK.RDK(r).rad_y RDK.RDK(r).rad_y],RDK.RDK(r).num,1)]'; % Generate initial dot position in a square space
    % check whether dots are within RDK area
    [dot_outside] = check_RDKarea(dot_pos,RDK.RDK(r));
    [dot_inside] = check_crosscutout(dot_pos,RDK.RDK(r),cross);
   
    
    while any([dot_outside.any; dot_inside]) % generates new positions for only the invalid dots
        dot_pos(:,dot_outside.any) = [randi([- RDK.RDK(r).rad_x  RDK.RDK(r).rad_x],sum(dot_outside.any),1) randi([- RDK.RDK(r).rad_y  RDK.RDK(r).rad_y],sum(dot_outside),1)]';
        dot_pos(:,dot_inside) = [randi([- RDK.RDK(r).rad_x  RDK.RDK(r).rad_x],sum(dot_inside),1) randi([- RDK.RDK(r).rad_y  RDK.RDK(r).rad_y],sum(dot_inside),1)]';
        [dot_outside] = check_RDKarea(dot_pos, RDK.RDK(r));
        [dot_inside] = check_crosscutout(dot_pos, RDK.RDK(r),cross);
    end
    
    % setup dot position matrix
    dotmat_r = zeros(2,  RDK.RDK(r).num, frames.pertrial); %setup dot matrix: xy x dot number x frames
    
    for f = 1:frames.pertrial %predefine dot positions for every frame
            
        %update dot position either random or as target event
        if any(f>=trial.event.onset & f<trial.event.onset+(RDK.event.duration*scr.refrate) & trial.event.RDK == r)
            % do coherent motion of certain dots
            t.evidx = find(f>=trial.event.onset & f<trial.event.onset+(RDK.event.duration*scr.refrate));
            if any(f==trial.event.onset(t.evidx))
                % index which dots are moved coherently
                t.cohidx = randsample(RDK.RDK(r).num,ceil(RDK.RDK(r).num*RDK.event.coherence));
                t.cohidx2 = false(1,RDK.RDK(r).num); t.cohidx2(t.cohidx)=true;              
            end
            % do event coherent movement
            dot_pos(:,t.cohidx2) = dot_pos(:,t.cohidx2) + RDK.event.direction(trial.event.direction(t.evidx),:)'*RDK.RDK(r).mov_speed*movdot(f);
            % do random movement for the rest
            % randomize directions first
            t.rdiridx = randi(length(RDK.RDK(r).mov_dir),1,sum(~t.cohidx2));
            % apply to position
            dot_pos(:,~t.cohidx2) = dot_pos(:,~t.cohidx2) + RDK.RDK(r).mov_dir(t.rdiridx,:)'*RDK.RDK(r).mov_speed*movdot(f);
            
        else
            % randomize directions first
            t.rdiridx = randi(length(RDK.RDK(r).mov_dir),1,RDK.RDK(r).num);
            % apply to position
            dot_pos = dot_pos + RDK.RDK(r).mov_dir(t.rdiridx,:)'*RDK.RDK(r).mov_speed*movdot(f);
        end
    
        % check whether dots are within RDK area
        [dot_outside,dot_crit] = check_RDKarea(dot_pos,RDK.RDK(r));
        % push outside dots to other side
        switch RDK.RDK(r).shape
            case 0 % ellipse
                dot_pos(:,dot_outside.ellipse) = fix(-dot_pos(:,dot_outside.ellipse).*(1+(1-dot_crit.xy(dot_outside.ellipse)))); % mirror dot to other side and push in a bit
            case 1 % square % case {1,2}??
                dot_pos(1,dot_outside.x) = fix(-dot_pos(1,dot_outside.x).*(1+(1-dot_crit.x(dot_outside.x)))); % for squares x dimension
                dot_pos(2,dot_outside.y) = fix(-dot_pos(2,dot_outside.y).*(1+(1-dot_crit.y(dot_outside.y)))); % for squares y dimension
            
        end

        % check whether dots are outside fixation cross area (if necessary)
        [dot_inside,cross_crit] = check_crosscutout(dot_pos,RDK.RDK(r),cross);
        % push outside dots to other side
        dot_pos(:,dot_inside) = fix(-dot_pos(:,dot_inside).*(1+(1-cross_crit(dot_inside))))+1; % dot mirrored to the other side of cross cutout and pushed out a bit

        
        
        %write into dot matrix
        dotmat_r(:,:,f) = dot_pos; 
    end
     
    % add QUAD4X shift to positions
    if Propixx == 4
        for p=1:Propixx
        shift_idx = p:Propixx:frames.pertrial;  
        dotmat_r(:,:,shift_idx) = [dotmat_r(1,:,shift_idx)+scr.shift(p,1); dotmat_r(2,:,shift_idx)+scr.shift(p,2)];
        end
    end   
    % add RDK centershift to positions
    dotmat_r = [dotmat_r(1,:,:)+RDK.RDK(r).centershift(1); dotmat_r(2,:,:)+RDK.RDK(r).centershift(2)];
     
    % setup array to mark to wich RDK a dot belongs %ToDo: is this necessary?
    rdkidx_r(1:RDK.RDK(r).num,1:frames.pertrial) = r;
    
    % define size for each dot
    dotsize_r = repmat(RDK.RDK(r).dot_size,size(rdkidx_r));
    
    % append everything; dots of different RDKs are merged into one dimension
    colmat = cat(2,colmat,colmat_r);
    dotmat = cat(2,dotmat,dotmat_r);
    rdkidx = cat(1,rdkidx,rdkidx_r);
    dotsize = cat(1,dotsize,dotsize_r);
    
end

%% Shuffle

%shuffle the dots of the different RDKs, so that the lastly drawn RDK is not systematically overlapping the others
shuffle_idx = Shuffle(1:size(rdkidx,1));
rdkidx = rdkidx(shuffle_idx,:); %ToDo: weg?
colmat = colmat(:,shuffle_idx,:);
dotmat = dotmat(:,shuffle_idx,:);
dotsize = dotsize(shuffle_idx,:);

colmat = reshape(colmat,[size(colmat,1), size(rdkidx,1)*Propixx, frames.flips]);
dotmat = reshape(dotmat,[2, size(rdkidx,1)*Propixx, frames.flips]);
rdkidx = reshape(rdkidx,[size(rdkidx,1)*Propixx, frames.flips]);
dotsize = reshape(dotsize,[size(dotsize,1)*Propixx, frames.flips]);


end

%% SUBFUNCTIONS

function [dot_outside,dot_crit] = check_RDKarea(dot_pos,RDK)
dot_crit.x = dot_pos(1,:).^2/RDK.rad_x^2;
dot_crit.y = dot_pos(2,:).^2/RDK.rad_y^2;
dot_crit.xy = dot_crit.x + dot_crit.y;

dot_outside.ellipse = dot_crit.xy > 1;
dot_outside.x = dot_crit.x > 1;
dot_outside.y = dot_crit.y > 1;
dot_outside.square = (dot_outside.x + dot_outside.y) >= 1;

switch RDK.shape
    case 1
    dot_outside.any = dot_outside.square;
    case 0
    dot_outside.any = dot_outside.ellipse;
end

end

function [dot_inside,cross_crit] = check_crosscutout(dot_pos,RDK,cross)
 cross_crit = (dot_pos(1,:)+RDK.centershift(1)).^2/cross.rad^2 ...
    + (dot_pos(2,:)+RDK.centershift(2)).^2/cross.rad^2; % dots inside the fixation cross cutout (if cross.cutout nonzero); Note: these calculations still assume that RDK is central, centershift is added later; that's why we need to add centershift here for checking;
dot_inside = cross_crit <= 1;
end
