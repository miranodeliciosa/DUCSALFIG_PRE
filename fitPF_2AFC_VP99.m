%% Psychometric Function fit (2AFC) with Palamedes: Weibull + MLE
% Recs used:
% - Fix gamma=0.5 (2AFC), free small lapse (0–0.05)
% - Goodness-of-fit (deviance) and parametric bootstrap CIs

clear; clc;

%% -------- Paths & Load data --------
% Palamedes in working directory (adjust if needed)
if exist('./Palamedes', 'dir') == 7
    addpath(genpath('./Palamedes'));
else
    warning('Palamedes folder not found in working directory. Adjust addpath below.');
    % addpath(genpath('/path/to/Palamedes'));  % <- adjust if needed
end

%dataFile = '/Users/sebastianwehle/Documents/MATLAB/Data_DuCSalFigPRE/logfiles/VP99_timing.mat';
dataFile = '/home/pc/matlab/user/sebastian/DuC_salientfigure/PreExp/logfiles/VP99_timing.mat';
load(dataFile);   % loads variable(s) incl. "resp"

%% -------- Aggregate behavioral data from ALL trials --------
allAnglesVals = [];          % physical angles per event (across all trials)
allIsCorrect  = [];          % correctness per event (1=correct, 0=incorrect)

for t = 1:numel(resp.experiment)
    E = resp.experiment{t};
    if isempty(E), continue; end

    % Optional override (if needed): inject mapping
    SBA.event.anglesToRotate = 10:5:45;   % your FOR CODING values
    E.SBA.event.anglesToRotate = SBA.event.anglesToRotate;

    % Map contrast index to physical angle
    if isfield(E, 'SBA') && isfield(E.SBA, 'event') && isfield(E.SBA.event, 'anglesToRotate')
        angles = E.SBA.event.anglesToRotate(:)';               % 1xN
    else
        angles = [5 8 10 12 14 16 20 35];                      % fallback
        warning('Trial %d: missing anglesToRotate, using fallback.', t);
    end

    % Get contrast indices and response class for this trial
    if ~isfield(E, 'SBA_contrast') || ~isfield(E, 'event_response_class')
        warning('Trial %d missing SBA_contrast or response data. Skipping.', t);
        continue;
    end

    contrastIdx = E.SBA_contrast(:);
    respClass   = string(E.event_response_class(:));

    % Basic alignment
    nEv = min(numel(contrastIdx), numel(respClass));
    contrastIdx = contrastIdx(1:nEv);
    respClass   = respClass(1:nEv);

    % Valid contrast indices
    validIdx = contrastIdx >= 1 & contrastIdx <= numel(angles);
    contrastIdx = contrastIdx(validIdx);
    respClass   = respClass(validIdx);

    % Map to angles
    anglesPerEvent = angles(contrastIdx);

    % Convert responses to binary
    isCorrect = NaN(numel(respClass), 1);
    isCorrect(contains(lower(respClass), 'correct'))   = 1;
    isCorrect(contains(lower(respClass), 'incorrect')) = 0;

    % Keep only valid responses
    valid = ~isnan(isCorrect);
    allAnglesVals = [allAnglesVals; anglesPerEvent(valid)'];
    allIsCorrect  = [allIsCorrect;  isCorrect(valid)];
end

%% -------- Aggregate per contrast level --------
[stimLevels, ~, grp] = unique(allAnglesVals);   % sorted unique angles
numLevels = numel(stimLevels);

outOfNum   = accumarray(grp, 1, [numLevels 1], @sum, 0);
numCorrect = accumarray(grp, allIsCorrect, [numLevels 1], @sum, 0);

% Remove empty levels (safety check)
haveData   = outOfNum > 0;
stimLevels = stimLevels(haveData);
outOfNum   = outOfNum(haveData);
numCorrect = numCorrect(haveData);

% Proportion correct
propCorrect = numCorrect ./ outOfNum;
% % % % % 
% % % % % E = resp.experiment{1};   % trial = 1 (as you indicated)
% % % % % %% -------- Extract per-event data --------
% % % % % % Required fields:
% % % % % %   - E.SBA_contrast: vector of contrast indices (1..8)
% % % % % %   - E.event_response_class: cell array of 'correct'/'incorrect'
% % % % % %   - E.SBA.event.anglesToRotate: [5 8 10 12 14 16 20 35] (mapping from index->angle)
% % % % % 
% % % % % % Ensure column vectors
% % % % % contrastIdx = E.SBA_contrast(:);                           % 1..8
% % % % % respClass   = string(E.event_response_class(:));           % "correct"/"incorrect"
% % % % % 
% % % % % % Map indices to physical contrast values (angles)
% % % % % if isfield(E, 'SBA') && isfield(E.SBA, 'event') && isfield(E.SBA.event, 'anglesToRotate')
% % % % %     angles = E.SBA.event.anglesToRotate(:)';               % 1x8
% % % % % else
% % % % %     % Fallback to the values you provided in the message:
% % % % %     angles = SBA.event.anglesToRotate(:);                      % degrees
% % % % %     warning('Using fallback anglesToRotate from SBA.event.anglesToRotate(:).');
% % % % % end
% % % % % 
% % % % % % Basic sanity checks
% % % % % assert(all(contrastIdx >= 1 & contrastIdx <= numel(angles)), 'Contrast indices out of range.');
% % % % % assert(numel(respClass) == numel(contrastIdx), 'Length mismatch: responses vs contrasts.');
% % % % % 
% % % % % % Convert responses to 1=correct, 0=incorrect; ignore anything else
% % % % % isCorrect = NaN(numel(respClass), 1);
% % % % % isCorrect(contains(lower(respClass), 'correct'))   = 1;
% % % % % isCorrect(contains(lower(respClass), 'incorrect')) = 0;
% % % % % 
% % % % % valid = ~isnan(isCorrect);
% % % % % contrastIdx = contrastIdx(valid);
% % % % % isCorrect   = isCorrect(valid);
% % % % % 
% % % % % %% -------- Aggregate per contrast level --------
% % % % % levels = 1:numel(angles);
% % % % % stimLevels = angles(levels);       % x-axis (physical units: degrees)
% % % % % 
% % % % % outOfNum  = zeros(size(levels));
% % % % % numCorrect = zeros(size(levels));
% % % % % 
% % % % % for k = levels
% % % % %     idx = contrastIdx == k;
% % % % %     outOfNum(k)  = sum(idx);
% % % % %     numCorrect(k) = sum(isCorrect(idx));
% % % % % end
% % % % % 
% % % % % % Remove empty levels (no trials), if any
% % % % % haveData = outOfNum > 0;
% % % % % stimLevels = stimLevels(haveData);
% % % % % outOfNum   = outOfNum(haveData);
% % % % % numCorrect = numCorrect(haveData);
% % % % % 
% % % % % % Proportion correct (for plotting)
% % % % % propCorrect = numCorrect ./ outOfNum;


%% -------- Palamedes setup: Weibull PF + MLE --------
PF = @PAL_Weibull;                % choices: PAL_Weibull, PAL_Logistic, PAL_CumulativeNormal, ...

% Parameters: [alpha beta gamma lambda]
%   alpha: threshold/location
%   beta: slope
%   gamma: guess rate (fixed 0.5 for 2AFC)
%   lambda: lapse rate
paramsFree = [1 1 0 1];           % free: alpha,beta,lambda; fix gamma
gammaFixed = 0.5;                 % 2AFC
lapseLimits = [0 0.05];           % allow small lapses

% Search grid (broad but sensible). Alpha near data range; beta log-spaced.
searchGrid.alpha  = linspace(min(stimLevels), max(stimLevels), 200);
searchGrid.beta   = logspace(-1, 2, 121);      % ~0.1 to 100
searchGrid.gamma  = gammaFixed;
searchGrid.lambda = linspace(lapseLimits(1), lapseLimits(2), 51);

% Fit by MLE
[paramsValues, LL, exitflag] = PAL_PFML_Fit( ...
    stimLevels, numCorrect, outOfNum, ...
    searchGrid, paramsFree, PF, ...
    'lapseLimits', lapseLimits, ...
    'gammaEQlambda', false);

alpha  = paramsValues(1);
beta   = paramsValues(2);
gamma  = paramsValues(3); % #ok<NASGU>  % should be 0.5
lambda = paramsValues(4);

fprintf('--- MLE Fit (Weibull) ---\n');
fprintf('alpha (threshold): %.4f\n', alpha);
fprintf('beta   (slope)   : %.4f\n', beta);
fprintf('gamma  (guess)   : %.4f (fixed)\n', searchGrid.gamma);
fprintf('lambda (lapse)   : %.4f\n', lambda);
fprintf('Log-likelihood   : %.4f\n', LL);
fprintf('Exitflag         : %d\n\n', exitflag);

%% -------- Goodness-of-fit (deviance) --------
% Monte Carlo deviance test (increase numMC for stricter test)
%% -------- Goodness-of-fit (deviance) --------
numMC = 1000;  % number of Monte-Carlo simulations

% [Dev, pDev] = PAL_PFML_GoodnessOfFit( ...
%     stimLevels, numCorrect, outOfNum, ...
%     paramsValues, paramsFree, PF, ...
%     'searchGrid', searchGrid, ...
%     'lapseLimits', lapseLimits, ...
%     'gammaEQlambda', false, ...
%     numMC);   % <-- B goes here, as the final positional arg
numMC = 1000;

[Dev, pDev] = PAL_PFML_GoodnessOfFit( ...
    stimLevels, numCorrect, outOfNum, ...
    paramsValues, paramsFree, numMC, PF, ...
    'searchGrid', searchGrid, ...
    'lapseLimits', lapseLimits, ...
    'gammaEQlambda', false, ...
    'noGrouping', true);  % <-- prevents GroupTrialsbyX from being called



fprintf('--- Goodness-of-Fit (Deviance) ---\n');
fprintf('Deviance: %.3f,  p-value: %.4f  (numMC = %d)\n\n', Dev, pDev, numMC);

%% -------- Parametric bootstrap for 95%% CIs --------
numSim = 1000;
[SD, paramsSim, LLSim, exitflagSim] = PAL_PFML_BootstrapParametric( ...
    stimLevels, outOfNum, paramsValues, paramsFree, numSim, PF, ...
    'lapseLimits', lapseLimits, ...
    'gammaEQlambda', false); %#ok<ASGLU,NASGU>

CI = prctile(paramsSim, [2.5 97.5], 1);  % rows: [2.5% 97.5%], cols: alpha beta gamma lambda
fprintf('--- 95%% Parametric Bootstrap CIs ---\n');
fprintf('alpha:  [%.4f, %.4f]\n', CI(1,1), CI(2,1));
fprintf('beta :  [%.4f, %.4f]\n', CI(1,2), CI(2,2));
fprintf('gamma:  fixed at 0.5\n');
fprintf('lambda: [%.4f, %.4f]\n\n', CI(1,4), CI(2,4));

%% -------- Plot: PF fit + empirical data --------
PFhat = @(x,a,b,g,l) PF([a b g l], x);
gammaFixed = searchGrid.gamma;  % 0.5 for 2AFC

% Generate smooth fit curve
xPlot = linspace(min(stimLevels), max(stimLevels), 400);
yFit  = PFhat(xPlot, alpha, beta, gammaFixed, lambda);

% Plot
figure('Color','w'); hold on;

% Empirical data: proportion correct with binomial SE
pc = propCorrect;
se = sqrt(pc .* (1 - pc) ./ outOfNum);
errorbar(stimLevels, pc, se, 'o', ...
    'MarkerFaceColor', [0.2 0.2 0.2], ...
    'LineStyle', 'none', ...
    'DisplayName', 'Data (\pm1 SE)');

% Fitted psychometric function
plot(xPlot, yFit, 'LineWidth', 2, 'DisplayName', 'Weibull fit');

% Aesthetics
ylim([0.45 1.01]);
xlim([min(stimLevels) max(stimLevels)]);
xlabel('Contrast (degrees)');
ylabel('Proportion correct');
title(sprintf('2AFC PF fit (Weibull, MLE)\n\\alpha = %.3f, \\beta = %.3f, \\lambda = %.3f', ...
    alpha, beta, lambda));
yline(0.5, ':', 'Chance', 'LabelHorizontalAlignment', 'left');
legend('Location', 'southeast');
grid on;

%% -------- Optional: threshold at 75% correct --------
crit = 0.75;
f = @(x) PFhat(x, alpha, beta, gammaFixed, lambda) - crit;

try
    x0 = alpha;  % starting point for fzero
    thr75 = fzero(f, x0);
    xline(thr75, '--', 'DisplayName', '75% threshold');
    text(thr75, 0.52, sprintf('%.2f° at 75%%', thr75), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
    fprintf('Threshold at %.0f%% correct: %.4f\n', 100 * crit, thr75);
catch
    warning('Threshold solver failed. Try checking the PF or expanding x-range.');
end
