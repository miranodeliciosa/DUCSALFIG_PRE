function F = MakeFlicker(Frames,CycleLength)

% Frames=0:(p.refreshrate*p.triallength)-1;
% CycleLength=p.refreshrate./p.frequencies;

% function F = MakeFlicker(Frames,CycleLength)
% Cycles start at Frame 0 (not Frame 1)
% "Frames" is a vector of all frame-numbers (e.g. 0:119)
% multiple CycleLength might be handed over (e.g. RefreshRate./[15 17.5])
% output F is then numel(CycleLength) x numel(Frames) with adjusted flicker
% intensities per frame
% In fact, this function realizes flicker interpolation according to Andersen & Müller (2015).
% Now quoting: "In the general case, our interpolation technique is
% defined as follows. With a monitor refresh rate R and a
% desired stimulus frequency f, the required cycle length ?
% is defined by: ? = R/f
% The stimulus intensity w (1:on, 0:off) in any given frame
% i can then be calculated as follows:
% w = { 1,                        if 1 <= i mod ? <= r_on*?
%       r_on*? + 1 ? i mod ?,     if r_on*? < i mod ? < r_on*? + 1
%       0,                        if r_on*? <= i mod ?
%       i mod ?,                  if i mod ? < 1
% where ron denotes the fraction of the stimulus cycle
% in which the stimulus is on (ron = 0.5 in the recordings
% reported here). Note that i modulo ? is the position
% within the current flicker cycle. The first line defines
% the frames in which the stimulus is on at full intensity,
% the second line the on–off transitions, the third line the
% off-frames and finally the fourth line defines the off–on
% transitions.  
% In order to present stimuli at intermediate intensities
% 0 < w <1 two further things need to be taken into account.
% First, the off-phase of a stimulus is not necessarily black
% but could be any arbitrary color C_off and thus intermediate
% stimulus intensities C must be a weighted average of
% the stimulus color C_on and the background or ‘off’ color
% C_off. Second, the relationship between color values in a
% stimulation program V_in (e.g. RGB values) and the output
% of a computer monitor V_out is not linear, but defined by a
% power function with an exponent gamma (?):
% V_out ~ V_in^?
% The exponent ? depends on the specific graphics hardware
% and its settings and usually lies between 1.8 and 2.2
% for standard computer equipment. Taking the two points
% above into account, the output color C can be computed
% as:
% C=( W * C1.^Gamma + (1-W) * C2.^Gamma ) .^(1/Gamma); "(end quoting)
% MakeFlicker will return (non)interpolated intensity values w (i.e. output F). 
% The actual color C have to be calculated elsewhere whith stimulus color
% and background color as further inputs, because these values will be
% linearly weighted according to w
%
% (c) and all credits to S. K. Andersen
% 2018 - NF commented on function usage
F=zeros(numel(CycleLength),numel(Frames));
for CycleNr=1:numel(CycleLength)
    FramesOn=CycleLength(CycleNr)/2; % "r_on*?"
    x=mod(Frames+1,CycleLength(CycleNr)); % "i mod ?"
    y=zeros(size(x)); % this will be the adjusted intensity vector w
    y(x<1)=x(x<1); % "w = i mod ?, if i mod ? < 1"
    y(x>=1 & x<=FramesOn)=1; % "w = 1, if 1 <= i mod ? <= r_on*?"
    y(x>FramesOn & x<FramesOn+1)=1+FramesOn-x(x>FramesOn & x<FramesOn+1); % "w = r_on*? + 1 ? i mod ?, if r_on*? < i mod ? < r_on*? + 1"
    y(x>=FramesOn+1)=0; % "w = 0, if r_on*? <= i mod ?"
    F(CycleNr,:)=y; % again w
end