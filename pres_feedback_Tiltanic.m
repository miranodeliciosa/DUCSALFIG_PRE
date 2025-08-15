function pres_feedback_DuCSalFig(responses,p,ps, key)
% PRES_FEEDBACK_DUCSALFIG  Displays subject performance feedback on screen.
%
%   PRES_FEEDBACK_DUCSALFIG(RESPONSES, P, PS, KEY, SBA, RDK)
%   computes the percentage of correct / incorrect / no-response (NoResp)
%   events and the overall mean reaction time (RT, in ms, across responded
%   events) from the RESPONSES structure. It prints these values to the
%   MATLAB console (in German) and presents them via Psychtoolbox on the
%   subject display until the key specified in KEY.SECRET is pressed.
%
%   SYNTAX
%     pres_feedback_DuCSalFig(responses, p, ps, key, SBA, RDK)
%
%   INPUTS
%     responses  Structure with the following fields:
%                - event_response_class : {nEvents x 1} cell array containing
%                    'correct' | 'incorrect' | [NaN] (NaN = no response)
%                - event_response_RT    : [nEvents x 1] double, reaction time
%                    in milliseconds; NaN for events without a response
%
%     p          (optionally used) structure; for text color:
%                - p.crs.color          : [1x3] uint8/double RGB (0–255).
%                  If missing, defaults to [220 220 220].
%
%     ps         Structure containing the already opened Psychtoolbox window:
%                - ps.window            : Window pointer from Screen('OpenWindow',...)
%
%     key        Structure with keyboard codes:
%                - key.SECRET           : Key code of the key that ends feedback
%                                         (e.g., KbName('q')).
%
%   OUTPUTS
%     (no return value) – The function produces:
%       - Console output with percentages and RT statistics (German text)
%       - On-screen feedback via Psychtoolbox until KEY.SECRET is pressed
%
%   CALCULATIONS
%     Percentages:
%       pctCorrect   = 100 * (# 'correct')   / nEvents
%       pctIncorrect = 100 * (# 'incorrect') / nEvents
%       pctNoResp    = 100 * (# NaN)         / nEvents
%
%     Reaction time:
%       rt_overall_mean = mean(event_response_RT, 'omitnan')
%       rt_overall_std  = std(event_response_RT, [], 'omitnan')
%     (NaN RT values are ignored; i.e., only responded events are averaged.)
%
%   SIDE EFFECTS
%     - Writes formatted feedback to the MATLAB console (fprintf).
%     - Draws text into the Psychtoolbox window PS.WINDOW and calls Screen('Flip',...).
%     - Blocks execution until KEY.SECRET is pressed (polling via KbCheck).
%
%   DEPENDENCIES
%     Requires Psychtoolbox (Screen, DrawFormattedText, KbCheck).
%     PS.WINDOW must already be opened via Screen('OpenWindow', ...).
%
%   ASSUMPTIONS / NOTES
%     - RESPONSES.EVENT_RESPONSE_CLASS can contain mixed cell content:
%       'correct' / 'incorrect' (char) and [NaN] (double) for NoResp.
%     - Times are in milliseconds.
%     - KEY.SECRET must be a valid index in KEYCODE (from KbCheck).
%
%   EXAMPLE
%     % Assuming ps.window is open, p.crs.color exists, and key.SECRET = KbName('q'):
%     pres_feedback_DuCSalFig(responses, p, ps, key, [], []);

%% average responses: 
%% ==== 1) Compute per-subject metrics from your new fields ====
% Supports either a single struct or an array: responses(1), responses(2), ...

isCorrect   = cellfun(@(x) ischar(x)    && strcmp(x,'correct'),   responses.event_response_class);
isIncorrect = cellfun(@(x) ischar(x)    && strcmp(x,'incorrect'), responses.event_response_class);
isNoResp    = cellfun(@(x) isnumeric(x) && isnan(x),               responses.event_response_class);

numEvents = numel(responses.event_response_class);
pctCorrect   = 100 * sum(isCorrect)   / numEvents;
pctIncorrect = 100 * sum(isIncorrect) / numEvents;
pctNoResp    = 100 * sum(isNoResp)    / numEvents;

% overall RT mean across all *responded* events
rt_overall_mean = mean(responses.event_response_RT, 'omitnan');   % use nanmean(rt) if needed
rt_overall_std = std(responses.event_response_RT, [], 'omitnan');   % use nanmean(rt) if needed

%% ==== 2) Console output (German) ====
fprintf(['\nKorrekt: %05.2f%% | Inkorrekt: %05.2f%% | Verpasst: %05.2f%%' ...
         ' | Reaktionszeit (M): %4.0f ms, (SD): %4.0f ms\n'], ...
        pctCorrect, pctIncorrect, pctNoResp, rt_overall_mean, rt_overall_std);
fprintf('\nMit "q" geht es weiter.\n');

%% ==== 3) On-screen (Psychtoolbox) feedback (German) ====
% We draw a single feedback page for one subject (sIdx). Change sIdx if you want a different subject.

% Choose text color (fall back to white if p/crs missing)
if exist('p','var') && isfield(p,'crs') && isfield(p.crs,'color')
    txtColor = p.crs.color;
else
    txtColor = [220 220 220];
end

% Build the German feedback text
text2present = [ ...
    'P A U S E' ...
    sprintf('\n\nValide: %1.0f%%',   pctCorrect+pctIncorrect) ...
    sprintf('\nKeine Antwort: %1.0f%%', pctNoResp) ...
    sprintf('\n\nReaktionszeit (M): %1.0f ms', rt_overall_mean) ...
];

% Draw
Screen('TextSize', ps.window, 36);
DrawFormattedText(ps.window, text2present, 'center', 'center', txtColor);

% Show it
Screen('Flip', ps.window);

% [key.pressed, key.firstPress]=KbQueueCheck;
key.rkey=key.SECRET;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    DrawFormattedText(ps.window, text2present, 'center', 'center', txtColor);
    Screen('Flip', ps.window, 0);                   % flip screen
end

end 