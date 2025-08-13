function event_onsets = generate_event_onsets_binned(nTrialsPerCondition, frameRate, flickerStep, nTimeBins, flag_plot)
% GENERATE_EVENT_ONSETS_BINNED 
% Erzeugt Event-Onsets mit Mindestabstand über gleichmäßig verteilte Zeitbins.
% Für Flicker-Frequenz 120/7 Hz (~17.14 Hz) → exakt 1 Onset alle 7 Frames.
%
% INPUT:
%   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–3)
%   frameRate           – z. B. 120
%   flickerStep         - z.B. 7 -> freq 120/7
%   nTimeBins           - Anzahl der Bins
%   flag_plot           – true/false für Visualisierung
%
% OUTPUT:
%   event_onsets – Zell-Array {1x3} mit je [nTrials x maxEvents] Eventframe-Matrizen

    % === Parameter ===
    numEvents = 4;
    %rng('default'); rng(randi());
    maxTries = 5000;
    minSeparation = 0.8;  % 800 ms
    minSeparationFrames = round(minSeparation * frameRate);
    binMargin = 0.1;

    % === Flicker-Frequenz 120/7 Hz → exakte Frame-Schritte ===
    flickerHz = frameRate / 7;
    %flickerStep = 7;  % 1 Onset alle 7 Frames
    totalDuration = 3.5;
    precuetime = 1.5;

    % === Erzeuge gleichmäßig verteilte mögliche Onset-Frames ===
    possibleFrames = 1:flickerStep:round(totalDuration * frameRate);
    totalFrames = numel(possibleFrames);
    frameTimes = (possibleFrames - 1) / frameRate - precuetime;

    % === Teile Frames in gleich große Gruppen (gleich viele pro Bin) ===
%     framesPerBin = floor(totalFrames / nTimeBins);
%     binFrameIdx = [1:framesPerBin:totalFrames];
%     binEdges = frameTimes(binFrameIdx);  % Zeitkanten basierend auf Frame-Zeit
%     binCenters = [binEdges(1:end-1) + diff(binEdges)/2 binEdges(end) + min(diff(binEdges)/2) ];
%     binIndexPerFrame = discretize(frameTimes, binEdges); % ?? NAN

    % === Teile Frames in gleich große Gruppen (gleich viele pro Bin) und zentriere um 0 ===
    framesPerBin = floor(totalFrames / nTimeBins);
    % 1. Finde den Frame, der Zeitpunkt 0 am nächsten liegt
    [~, zeroIdx] = min(abs(frameTimes));  % Index in possibleFrames
    % 2. Setze diesen als Bin-Mitte (oder -Start) und berechne Bin-Indizes symmetrisch
    halfBins = floor(nTimeBins / 2);
    startIdx = zeroIdx - framesPerBin * halfBins;
    binFrameIdx = startIdx + (0:framesPerBin:(framesPerBin * nTimeBins));
    % 3. Begrenzung auf gültige Indizes
    binFrameIdx = max(1, min(binFrameIdx, totalFrames));
    binFrameIdx = unique([binFrameIdx, totalFrames]);  % sicherstellen, dass letzte Grenze dabei ist
    % 4. Bin-Zeitkanten
    binEdges = frameTimes(binFrameIdx);
    binCenters = binEdges(1:end-1) + diff(binEdges)/2;
    % 5. Frame → Bin Mapping
    binIndexPerFrame = discretize(frameTimes, binEdges);


    % === Zielverteilung ===
    totalEvents = nTrialsPerCondition * sum(1:numEvents);
    idealPerBin = round(totalEvents / nTimeBins);
    maxPerBin = ceil(idealPerBin * (1 + binMargin));
    minPerBin = floor(idealPerBin * (1 - binMargin));
    binCounts = zeros(1, nTimeBins);

    % === Priorität für 1-Event-Trials (optional) ===
    priorityWeights = normpdf(binCenters, 0.3, 0.25);  % Fokus auf 300ms
    priorityWeights = priorityWeights / max(priorityWeights);

    % === Generiere Trials ===
    event_onsets = cell(1, numEvents);
    usedFrames = false(1, totalFrames);
    allEventFrames = [];

    for nEvents = numEvents:-1:1
        trials = zeros(nTrialsPerCondition, numEvents);
        for t = 1:nTrialsPerCondition
            success = false;
            tries = 0;
            while ~success && tries < maxTries
                tries = tries + 1;

                rawWeights = max(idealPerBin - binCounts', 0);
                if nEvents == 1
                    rawWeights = rawWeights .* priorityWeights';
                end

                weights = zeros(1, totalFrames);
                validIdx = ~isnan(binIndexPerFrame);
                weights(validIdx) = rawWeights(binIndexPerFrame(validIdx));
                weights(usedFrames) = 0;
                if all(weights == 0)
                    warning('Keine gültigen Frames mehr verfügbar.');
                    break;
                end
                weights = weights / sum(weights);
                eligibleFrames = find(weights > 0);
                eligibleWeights = weights(eligibleFrames);
                eligibleWeights = eligibleWeights / sum(eligibleWeights);

                chosenIdx = randsample(length(eligibleFrames), nEvents, true, eligibleWeights);
                chosenFrames = possibleFrames(eligibleFrames(chosenIdx));
                candidate = sort(unique(chosenFrames));

                if all(diff(candidate) >= minSeparationFrames) && length(candidate) == nEvents
                    trials(t, 1:nEvents) = candidate;
                    usedFrames(candidate) = true;
                    allEventFrames = [allEventFrames, candidate];

                    % Bin-Zuordnung aktualisieren
                    candidate_idx = zeros(size(candidate));
                    for i = 1:numel(candidate)
                        candidate_idx(i) = find(possibleFrames == candidate(i), 1);
                    end
                    binsUsed = binIndexPerFrame(candidate_idx);
                    for b = binsUsed
                        if ~isnan(b)
                            binCounts(b) = binCounts(b) + 1;
                        end
                    end

                    success = true;
                end
            end
        end
        event_onsets{nEvents} = trials;
    end


    % === Visualisierung ===
    if flag_plot>0
        figure(flag_plot);
        colors = lines(numEvents);
        allTimes = [];

        for nEvents = 1:numEvents
            subplot(3, 2, nEvents);
            onsetFrames = event_onsets{nEvents};
            onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
            for i = 1:nTrialsPerCondition
                x = onsetTimesSec(i, 1:nEvents);
                y = i * ones(1, nEvents);
                scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
                allTimes = [allTimes, x];
            end
            title(sprintf('%d Event(s) pro Trial', nEvents));
            xlabel('Zeit (s)');
            ylabel('Trial');
            xlim([min(frameTimes), max(frameTimes)]);
            ylim([0 nTrialsPerCondition + 1]);
            grid on;
        end

        subplot(3, 2, [5]);
        histogram(allTimes, 'BinEdges', binEdges, ...
                  'FaceColor', [0.2 0.4 0.8], 'Normalization', 'count');
        hold on;
        plot(binCenters, idealPerBin, '--r');
        plot(binCenters, maxPerBin, ':r');
        plot(binCenters, minPerBin, ':r');
        xlabel('Zeit (s)');
        ylabel('Eventanzahl pro Bin');
        title('Globale Verteilung über Zeitbins');
        subtitle(sprintf('Binbreite: %.0f ms (variabel in Zeit)', 1000 * max(diff(binEdges))));
        xlim([min(frameTimes), max(frameTimes)]);
        grid on;

      
        subplot(3, 2, [6]);
        bar(binCenters, priorityWeights, 0.9, 'FaceColor', [0.3 0.6 0.3]);
        xlabel('Zeit (s)');
        ylabel('Priority Weight');
        title('Prioritätsgewichtung für 1-Event-Trials');
        grid on;
        xlim([min(frameTimes), max(frameTimes)]);

         sgtitle(sprintf('Events mit Mindestabstand (%.0f ms) bei %.2f Hz Flicker', ...
            1000 * minSeparation, flickerHz));

    % === Bin-Statistik ===
    tooLowBins = find(binCounts < minPerBin');
    tooHighBins = find(binCounts > maxPerBin');
    if ~isempty(tooLowBins)
        warning('⚠️  %d Bins sind unterbesetzt (< min): %s', ...
                length(tooLowBins), mat2str(tooLowBins));
    end
    if ~isempty(tooHighBins)
        warning('⚠️  %d Bins sind überbesetzt (> max): %s', ...
                length(tooHighBins), mat2str(tooHighBins));
    end

    fprintf('\n=== Bin-Auslastung (Frames gleich verteilt) ===\n');
    disp(table((1:nTimeBins)', binCounts', repmat(idealPerBin,nTimeBins,1) , repmat(minPerBin,nTimeBins,1), repmat(maxPerBin,nTimeBins,1), ...
        'VariableNames', {'Bin', 'EventCount', 'Ideal', 'MinTol', 'MaxTol'}));
    end
   
end

% function event_onsets = generate_event_onsets_binned(nTrialsPerCondition, frameRate, nTimeBins, flag_plot)
% % GENERATE_EVENT_ONSETS_BINNED 
% % Events werden über Zeitbins möglichst gleich verteilt, mit 800ms Mindestabstand.
% % 1-Event-Trials werden bevorzugt in einen priorisierten Zeitbereich gelegt.
% %
% % INPUT:
% %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–3)
% %   frameRate           – Hz (z. B. 120)
% %   nTimeBins           – Anzahl gewünschter Zeitbins
% %   flag_plot           – true/false für Visualisierung
% %
% % OUTPUT:
% %   event_onsets – Zell-Array {1x3} mit je [nTrials x maxEvents] Eventframe-Matrizen
% 
%     numEvents = 4;
%     rng('default');
%     rng(2);
%     maxTries = 5000;
%     binMargin = 0.1;
% 
%     freq = 120/7;
%     precuetime = 1.7;
%     totalDuration = 3.2;
%     minSeparationFrames = round(0.8 * frameRate);
% 
%     % === Nur Flicker-Onset-Frames zulassen ===
%     possibleFrames = ceil(1 : frameRate/freq : totalDuration * frameRate);
%     totalFrames = numel(possibleFrames);
%     timeVec = (possibleFrames - 1) / frameRate - precuetime;
% 
%     % === Zeitbins vorbereiten ===
%     binWidth = totalDuration / nTimeBins;
%     binStart = -precuetime;
%     binEnd = binStart + totalDuration;
%     offsetToZero = mod(-binStart, binWidth);
%     binEdges = (binStart + offsetToZero) : binWidth : binEnd;
%     if length(binEdges) < nTimeBins + 1
%         binEdges(end+1) = binEdges(end) + binWidth;
%     end
%     binIndexPerFrame = discretize(timeVec, binEdges);
%     binCenters = binEdges(1:end-1) + diff(binEdges)/2;
% 
%     % === Zielverteilung und Prioritätsgewichtung ===
%     totalEvents = nTrialsPerCondition * sum(1:numEvents);
%     availableFramesPerBin = accumarray(binIndexPerFrame(~isnan(binIndexPerFrame))', ...
%                                        1, [nTimeBins, 1]);
%     proportionalWeights = availableFramesPerBin / sum(availableFramesPerBin);
%     idealPerBin = round(totalEvents * proportionalWeights);
%     maxPerBin = ceil(idealPerBin .* (1 + binMargin));
%     minPerBin = floor(idealPerBin .* (1 - binMargin));
%     binCounts = zeros(1, nTimeBins);
% 
%     % Prioritätsverteilung für 1-Event-Trials (z. B. Peak bei 300 ms)
%     priorityWeights = normpdf(binCenters, 0.3, 0.25);
%     priorityWeights = priorityWeights / max(priorityWeights);
% 
%     % === Trialweise generieren ===
%     event_onsets = cell(1, numEvents);
%     usedFrames = false(1, totalFrames);
%     allEventFrames = [];
% 
%     for nEvents = numEvents:-1:1
%         trials = zeros(nTrialsPerCondition, numEvents);
%         for t = 1:nTrialsPerCondition
%             success = false;
%             tries = 0;
%             while ~success && tries < maxTries
%                 tries = tries + 1;
% 
%                 rawWeights = max(idealPerBin - binCounts', 0);
%                 if nEvents == 1
%                     rawWeights = rawWeights .* priorityWeights';
%                 end
% 
%                 weights = zeros(1, totalFrames);
%                 validIdx = ~isnan(binIndexPerFrame);
%                 weights(validIdx) = rawWeights(binIndexPerFrame(validIdx));
%                 weights(usedFrames) = 0;
%                 if all(weights == 0)
%                     warning('Keine gültigen Frames mehr verfügbar.');
%                     break;
%                 end
%                 weights = weights / sum(weights);
%                 eligibleFrames = find(weights > 0);
%                 eligibleWeights = weights(eligibleFrames);
%                 eligibleWeights = eligibleWeights / sum(eligibleWeights);
% 
%                 chosenIdx = randsample(length(eligibleFrames), nEvents, true, eligibleWeights);
%                 chosenFrames = possibleFrames(eligibleFrames(chosenIdx));
%                 candidate = sort(unique(chosenFrames));
% 
%                 if all(diff(candidate) >= minSeparationFrames) && length(candidate) == nEvents
%                     trials(t, 1:nEvents) = candidate;
%                     usedFrames(candidate) = true;
%                     allEventFrames = [allEventFrames, candidate];
% 
%                     % Bin-Zuordnung
%                     candidate_idx = zeros(size(candidate));
%                     for ican = 1:numel(candidate)
%                         candidate_idx(ican) = find(possibleFrames == candidate(ican), 1);
%                     end
%                     binsUsed = binIndexPerFrame(candidate_idx);
%                     for b = binsUsed
%                         if ~isnan(b)
%                             binCounts(b) = binCounts(b) + 1;
%                         end
%                     end
% 
%                     success = true;
%                 end
%             end
% 
%             if ~success
%                 warning('Trial %d mit %d Events konnte nach %d Versuchen nicht erstellt werden.', ...
%                         t, nEvents, maxTries);
%             end
%         end
%         event_onsets{nEvents} = trials;
%     end
% 
%     % === Visualisierung ===
%     if flag_plot
%         figure;
%         colors = lines(numEvents);
%         allTimes = [];
% 
%         for nEvents = 1:numEvents
%             subplot(3, 2, nEvents);
%             onsetFrames = event_onsets{nEvents};
%             onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
%             for i = 1:nTrialsPerCondition
%                 x = onsetTimesSec(i, 1:nEvents);
%                 y = i * ones(1, nEvents);
%                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
%                 allTimes = [allTimes, x];
%             end
%             title(sprintf('%d Event(s) pro Trial', nEvents));
%             xlabel('Zeit (s)');
%             ylabel('Trial');
%             xlim([min(timeVec), max(timeVec)]);
%             ylim([0 nTrialsPerCondition + 1]);
%             grid on;
%         end
% 
%         % Globale Histogramm-Verteilung
%         subplot(3, 2, [5 6]);
%         histogram(allTimes, 'BinEdges', binEdges, ...
%                   'FaceColor', [0.2 0.4 0.8], 'Normalization', 'count');
%         hold on;
%         plot(binCenters, idealPerBin, '--r');
%         plot(binCenters, maxPerBin, ':r');
%         plot(binCenters, minPerBin, ':r');
%         xlabel('Zeit (s)');
%         ylabel('Eventanzahl pro Bin');
%         title('Globale Verteilung über Zeitbins');
%         subtitle(sprintf('Binbreite: %.0f ms', 1000 * binWidth));
%         xlim([min(timeVec), max(timeVec)]);
%         grid on;
% 
%         sgtitle(sprintf('Events mit Mindestabstand (800 ms) über %d Zeitbins', nTimeBins));
% 
%         % Zusatzplot: Priorität
%         figure;
%         bar(binCenters, priorityWeights, 0.9, 'FaceColor', [0.3 0.6 0.3]);
%         xlabel('Zeit (s)');
%         ylabel('Priorität (1-Event-Trials)');
%         title('Priorisierte Bin-Zielverteilung (nur für 1-Event-Trials)');
%         grid on;
%     end
% 
%     % === Bin-Statistik ===
%     tooLowBins = find(binCounts < minPerBin');
%     tooHighBins = find(binCounts > maxPerBin');
%     if ~isempty(tooLowBins)
%         warning('⚠️  %d Bins sind unterbesetzt (< min): %s', ...
%                 length(tooLowBins), mat2str(tooLowBins));
%     end
%     if ~isempty(tooHighBins)
%         warning('⚠️  %d Bins sind überbesetzt (> max): %s', ...
%                 length(tooHighBins), mat2str(tooHighBins));
%     end
% 
%     fprintf('\n=== Bin-Auslastung (angepasst & priorisiert) ===\n');
%     disp(table((1:nTimeBins)', binCounts', idealPerBin, minPerBin, maxPerBin, ...
%         'VariableNames', {'Bin', 'EventCount', 'Ideal', 'MinTol', 'MaxTol'}));
% end
% 
% 
% % function [event_onsets, binstats] = generate_event_onsets_binned(nTrialsPerCondition, frameRate, nTimeBins, flag_plot)
% % % GENERATE_EVENT_ONSETS_BINNED 
% % % Events werden über Zeitbins möglichst gleich verteilt, mit 800ms Mindestabstand.
% % %
% % % INPUT:
% % %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–4)
% % %   frameRate           – Hz (z. B. 120)
% % %   nTimeBins           – Anzahl gewünschter Zeitbins
% % %   flag_plot           – true/false für Visualisierung
% % %
% % % OUTPUT:
% % %   event_onsets – Zell-Array {1x4} mit je [nTrials x maxEvents] Eventframe-Matrizen
% % 
% %     numEvents = 4;
% %     rng('default');
% %     rng(2);
% %     maxTries = 5000;
% %     binMargin = 0.1;
% % 
% %     freq = 17.8;
% %     precuetime = 1.7;
% %     totalDuration = 3.2;
% %     minSeparationFrames = round(0.8 * frameRate);
% % 
% %     % === Nur Flicker-Onset-Frames zulassen ===
% %     possibleFrames = ceil(1 : frameRate/freq : totalDuration * frameRate);
% %     totalFrames = numel(possibleFrames);
% %     timeVec = (possibleFrames - 1) / frameRate - precuetime;
% % 
% %     % === Zeitbins vorbereiten (auf 0 s zentriert) ===
% %     binWidth = totalDuration / nTimeBins;
% %     binStart = -precuetime;
% %     binEnd = binStart + totalDuration;
% %     offsetToZero = mod(-binStart, binWidth);
% %     binEdges = (binStart + offsetToZero) : binWidth : binEnd;
% %     if length(binEdges) < nTimeBins + 1
% %         binEdges(end+1) = binEdges(end) + binWidth;
% %     end
% %     binIndexPerFrame = discretize(timeVec, binEdges);
% % 
% %     % === Adaptive Zielverteilung ===
% %     totalEvents = nTrialsPerCondition * sum(1:numEvents);
% %     availableFramesPerBin = accumarray(binIndexPerFrame(~isnan(binIndexPerFrame))', ...
% %                                        1, [nTimeBins, 1]);
% %     proportionalWeights = availableFramesPerBin / sum(availableFramesPerBin);
% %     idealPerBin = round(totalEvents * proportionalWeights);
% %     maxPerBin = ceil(idealPerBin .* (1 + binMargin));
% %     minPerBin = floor(idealPerBin .* (1 - binMargin));
% %     binCounts = zeros(1, nTimeBins);
% % 
% %     % === Trialweise generieren ===
% %     event_onsets = cell(1, numEvents);
% %     usedFrames = false(1, totalFrames);
% %     allEventFrames = [];
% % 
% %     for nEvents = numEvents:-1:1
% %         trials = zeros(nTrialsPerCondition, numEvents);
% %         for t = 1:nTrialsPerCondition
% %             success = false;
% %             tries = 0;
% %             while ~success && tries < maxTries
% %                 tries = tries + 1;
% % 
% %                 rawWeights = max(idealPerBin - binCounts', 0);
% %                 weights = zeros(1, totalFrames);
% %                 validIdx = ~isnan(binIndexPerFrame);
% %                 weights(validIdx) = rawWeights(binIndexPerFrame(validIdx));
% %                 weights(usedFrames) = 0;
% %                 if all(weights == 0)
% %                     warning('Keine gültigen Frames mehr verfügbar.');
% %                     break;
% %                 end
% %                 weights = weights / sum(weights);
% %                 eligibleFrames = find(weights > 0);
% %                 eligibleWeights = weights(eligibleFrames);
% %                 eligibleWeights = eligibleWeights / sum(eligibleWeights);
% % 
% %                 try 
% %                     candidate = sort(unique(randsample(possibleFrames(eligibleFrames), nEvents, true, eligibleWeights)));
% %                 catch ME
% %                     ME.message
% %                     nEvents
% %                     t
% %                     candidate
% %                     warning('Nehme zufällige Zuordnung ohne Gewichtung für diesen Trial.')
% %                     eligibleFrames
% %                     candidate = sort(unique(randsample(possibleFrames, nEvents)));
% %                     candidate
% %                 end 
% % 
% %                 if all(diff(candidate) >= minSeparationFrames) && length(unique(candidate)) == nEvents
% %                     trials(t, 1:nEvents) = candidate;
% %                     usedFrames(candidate) = true;
% %                     allEventFrames = [allEventFrames, candidate];
% % 
% %                     % Bin-Zuordnung
% %                     candidate_idx = zeros(size(candidate));
% %                     for ican = 1:numel(candidate)
% %                         candidate_idx(ican) = find(possibleFrames == candidate(ican), 1);
% %                     end
% %                     binsUsed = binIndexPerFrame(candidate_idx);
% %                     for b = binsUsed
% %                         if ~isnan(b)
% %                             binCounts(b) = binCounts(b) + 1;
% %                         end
% %                     end
% % 
% %                     success = true;
% %                 end
% %             end
% % 
% %             if ~success
% %                 warning('Trial %d mit %d Events konnte nach %d Versuchen nicht erstellt werden. Zufällige Zurordnung .', ...
% %                         t, nEvents, maxTries);
% %                 candidate = sort(unique(randsample(possibleFrames(eligibleFrames), nEvents)));
% %                 if all(diff(candidate) >= minSeparationFrames) && length(unique(candidate)) == nEvents
% %                     trials(t, 1:nEvents) = candidate;
% %                     usedFrames(candidate) = true;
% %                     allEventFrames = [allEventFrames, candidate];
% %                 end
% %             end
% %         end
% %         event_onsets{nEvents} = trials;
% %     end
% % 
% %     % === Visualisierung ===
% %     if flag_plot
% %         figure;
% %         colors = lines(numEvents);
% %         allTimes = [];
% % 
% %         for nEvents = 1:numEvents
% %             subplot(3, 2, nEvents);
% %             onsetFrames = event_onsets{nEvents};
% %             onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
% %             for i = 1:nTrialsPerCondition
% %                 x = onsetTimesSec(i, 1:nEvents);
% %                 y = i * ones(1, nEvents);
% %                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
% %                 allTimes = [allTimes, x];
% %             end
% %             title(sprintf('%d Event(s) pro Trial', nEvents));
% %             xlabel('Zeit (s)');
% %             ylabel('Trial');
% %             xlim([min(timeVec), max(timeVec)]);
% %             ylim([0 nTrialsPerCondition + 1]);
% %             grid on;
% %         end
% % 
% %         subplot(3, 2, [5 6]);
% %         histogram(allTimes, 'BinEdges', binEdges, ...
% %                   'FaceColor', [0.2 0.4 0.8], 'Normalization', 'count');
% %         hold on;
% %         plot(binEdges(1:end-1) + binWidth/2, idealPerBin, '--r');
% %         plot(binEdges(1:end-1) + binWidth/2, maxPerBin, ':r');
% %         plot(binEdges(1:end-1) + binWidth/2, minPerBin, ':r');
% %         xlabel('Zeit (s)');
% %         ylabel('Eventanzahl pro Bin');
% %         title('Globale Verteilung über Zeitbins');
% %         subtitle(sprintf('Binbreite: %.0f ms', 1000 * binWidth));
% %         xlim([min(timeVec), max(timeVec)]);
% %         grid on;
% % 
% %         sgtitle(sprintf('Events mit Mindestabstand (800 ms) über %d Zeitbins', nTimeBins));
% %     end
% % 
% %     % === Bin-Statistik ===
% %     tooLowBins = find(binCounts < minPerBin');
% %     tooHighBins = find(binCounts > maxPerBin');
% %     if ~isempty(tooLowBins)
% %         warning('⚠️  %d Bins sind unterbesetzt (< min): %s', ...
% %                 length(tooLowBins), mat2str(tooLowBins));
% %     end
% %     if ~isempty(tooHighBins)
% %         warning('⚠️  %d Bins sind überbesetzt (> max): %s', ...
% %                 length(tooHighBins), mat2str(tooHighBins));
% %     end
% % 
% %     fprintf('\n=== Bin-Auslastung (angepasst an Frameverfügbarkeit) ===\n');
% %     binstats = table((1:nTimeBins)', binCounts', idealPerBin, minPerBin, maxPerBin,...
% %         'VariableNames', {'Bin', 'EventCount', 'Ideal', 'MinTol', 'MaxTol'});
% %     disp(binstats);
% % end
% %  
% % 
% % % function event_onsets = generate_event_onsets_binned(nTrialsPerCondition, frameRate, nTimeBins, flag_plot)
% % % % GENERATE_EVENT_ONSETS_BINNED 
% % %     %!!!! Funktioniert für onsets_binned = generate_event_onsets_binned(25, 120, 15, 1);
% % % % Events werden über Zeitbins möglichst gleich verteilt, mit 800ms Mindestabstand.
% % % %
% % % % INPUT:
% % % %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–4)
% % % %   frameRate           – Hz (z. B. 120)
% % % %   nTimeBins           – Anzahl gewünschter Zeitbins
% % % %   flag_plot           – true/false für Visualisierung
% % % %
% % % % OUTPUT:
% % % %   event_onsets – Zell-Array {1x4} mit je [nTrials x maxEvents] Eventframe-Matrizen
% % % 
% % % 
% % % 
% % % 
% % % 
% % %     numEvents = 3;
% % % 
% % %     rng('default');
% % %     rng(2);
% % %     maxTries = 5000;
% % %     binMargin = 0.1;  % ±5% Toleranz bei Zielanzahl pro Bin
% % % 
% % %     freq = 17.8;
% % %     framerate = 120;
% % %     precuetime = 1.7;
% % %     totalDuration = 3.2;
% % % %    totalFrames = round(totalDuration * frameRate);
% % %     minSeparationFrames = round(0.8 * frameRate);
% % %     %possibleFrames = 1:totalFrames;
% % %     possibleFrames = ceil(1:framerate/freq:totalDuration*framerate);
% % %     totalFrames = numel(possibleFrames);
% % %     timeVec = (possibleFrames - 1) / frameRate - precuetime;
% % % 
% % %     totalEvents = nTrialsPerCondition * sum(1:numEvents);
% % % 
% % %     % === Zeitbins vorbereiten ===
% % %         % bins auf Null zentrieren
% % %         binWidth = totalDuration / nTimeBins;
% % %         binStart = -precuetime;
% % %         binEnd   = binStart + totalDuration;
% % %         offsetToZero = mod(-binStart, binWidth);
% % %         binEdges = (binStart + offsetToZero) : binWidth : binEnd;
% % %         
% % %         % Fallback: falls Rundungsfehler einen Bin zu wenig macht
% % %         if length(binEdges) < nTimeBins + 1
% % %             binEdges(end+1) = binEdges(end) + binWidth;
% % %         end
% % % 
% % %     binCounts = zeros(1, nTimeBins);
% % %     binIndexPerFrame = discretize(timeVec, binEdges);
% % % 
% % %     framesPerBin = accumarray(binIndexPerFrame(~isnan(binIndexPerFrame))', 1, [nTimeBins, 1]);
% % %     totalAvailable = sum(framesPerBin);
% % %     idealPerBin = round(totalEvents * (framesPerBin / totalAvailable));
% % %     maxPerBin = ceil(idealPerBin * (1 + binMargin));
% % %     minPerBin = floor(idealPerBin * (1 - binMargin));
% % % 
% % % 
% % % 
% % % % DEBUGGING
% % % % counts = accumarray(binIndexPerFrame(~isnan(binIndexPerFrame))', 1, [nTimeBins, 1]);
% % % % figure; 
% % % % bar(counts);
% % % % title('Verfügbarkeit von Flicker-Onset-Frames pro Zeitbin');
% % % % xlabel('Bin');
% % % % ylabel('Frame Count');
% % % 
% % % 
% % %     % === Trialweise generieren ===
% % % 
% % %     event_onsets = cell(1, numEvents);
% % %     usedFrames = false(1, totalFrames);
% % %     allEventFrames = [];
% % % 
% % % 
% % % 
% % % 
% % %     for nEvents = numEvents:-1:1
% % %         trials = zeros(nTrialsPerCondition, numEvents);
% % %         for t = 1:nTrialsPerCondition
% % %             success = false;
% % %             tries = 0;
% % %             while ~success && tries < maxTries
% % %                 tries = tries + 1;
% % % 
% % %                 % Gewichtete Auswahl basierend auf Unterbesetzung
% % %                 rawWeights = max(idealPerBin - binCounts, 0);
% % %                 weights = zeros(1, totalFrames);
% % %                 validIdx = ~isnan(binIndexPerFrame);
% % %                 weights(validIdx) = rawWeights(binIndexPerFrame(validIdx));
% % %                 weights(usedFrames) = 0;
% % %                 if all(weights == 0)
% % %                     warning('Keine gültigen Frames mehr verfügbar. Weise Frame zufällig zu.');
% % % 
% % %                     break;
% % %                 end
% % %                 weights = weights / sum(weights);
% % %                 eligibleFrames = find(weights > 0);
% % %                 % Nur gültige Frames herausfiltern
% % %                 eligibleWeights = weights(eligibleFrames);  % Jetzt hat gleiche Länge wie eligibleFrames
% % %                 eligibleWeights = eligibleWeights / sum(eligibleWeights);  % Normieren (sicherheitsbedingt)
% % % 
% % %                 % Dann samplen
% % %                 candidate = sort(unique(randsample(possibleFrames(eligibleFrames), nEvents, true, eligibleWeights)));
% % % 
% % % 
% % %                 % Prüfe Mindestabstand
% % %                 if all(diff(candidate) >= minSeparationFrames) && length(unique(candidate)) == nEvents
% % %                     trials(t, 1:nEvents) = candidate;
% % %                     usedFrames(candidate) = true;
% % %                     allEventFrames = [allEventFrames, candidate];
% % %                 
% % %                     % Aktualisiere Bin-Zähler
% % %                     for ican= 1:numel(candidate)
% % %                         candidate_idx(ican) = find(possibleFrames(eligibleFrames) == candidate(ican));
% % %                     end 
% % %                     binsUsed = binIndexPerFrame(candidate_idx);
% % %                     for b = binsUsed
% % %                         if ~isnan(b)
% % %                             binCounts(b) = binCounts(b) + 1;
% % %                         end
% % %                     end
% % %                 
% % %                     success = true;
% % %                 end
% % %             end
% % % 
% % %             if ~success
% % %                 warning('Trial %d mit %d Events konnte nach %d Versuchen nicht erstellt werden. Letzte Kandidaten werden trotzdem verwendet.', ...
% % %                         t, nEvents, maxTries);
% % % %                 % Alternativ: Die Bedingung der Gleichverteilung
% % % %                 entschärfen und zufällig ziehen 
% % %                   candidate = sort(unique(randsample(possibleFrames(eligibleFrames), nEvents)));
% % %                     if all(diff(candidate) >= minSeparationFrames) && length(unique(candidate)) == nEvents
% % %                     trials(t, 1:nEvents) = candidate;
% % %                     usedFrames(candidate) = true;
% % %                     allEventFrames = [allEventFrames, candidate];
% % %                     end 
% % % %                 trials(t, 1:nEvents) = candidate; 
% % % %                 usedFrames(candidate) = true;
% % % %                 allEventFrames = [allEventFrames, candidate];
% % % %                 binsUsed = binIndexPerFrame(candidate);
% % % %                 for b = binsUsed
% % % %                     if ~isnan(b)
% % % %                         binCounts(b) = binCounts(b) + 1;
% % % %                     end
% % % %                 end
% % %             end
% % %         end
% % %         event_onsets{nEvents} = trials;
% % %     end
% % % 
% % %     % === Visualisierung ===
% % %     if flag_plot
% % %         figure;
% % %         colors = lines(numEvents);
% % %         allTimes = [];
% % % 
% % %         for nEvents = 1:numEvents
% % %             subplot(3, 2, nEvents);
% % %             onsetFrames = event_onsets{nEvents};
% % %             onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
% % %             for i = 1:nTrialsPerCondition
% % %                 x = onsetTimesSec(i, 1:nEvents);
% % %                 y = i * ones(1, nEvents);
% % %                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
% % %                 allTimes = [allTimes, x];
% % %             end
% % %             title(sprintf('%d Event(s) pro Trial', nEvents));
% % %             xlabel('Zeit (s)');
% % %             ylabel('Trial');
% % %             xlim([min(timeVec), max(timeVec)]);
% % %             ylim([0 nTrialsPerCondition + 1]);
% % %             grid on;
% % %         end
% % % 
% % %         subplot(3, 2, [5 6]);
% % %         histogram(allTimes, 'BinEdges', binEdges, ...
% % %                   'FaceColor', [0.2 0.4 0.8], 'Normalization', 'count');
% % %         hold on;
% % %         yline(idealPerBin, '--r', 'Ziel');
% % %         yline(maxPerBin, ':r', 'Max Toleranz');
% % %         yline(minPerBin, ':r', 'Min Toleranz');
% % %         xlabel('Zeit (s)');
% % %         ylabel('Eventanzahl pro Bin');
% % %         title('Globale Verteilung über Zeitbins');
% % %         subtitle(sprintf('Binbreite: %.0f ms', 1000 * min(diff(binEdges))));
% % %         xlim([min(timeVec), max(timeVec)]);
% % %         grid on;
% % % 
% % %         sgtitle(sprintf('Events mit Mindestabstand (800 ms) über %d Zeitbins', nTimeBins));
% % %     end
% % % 
% % %     % === Bin-Statistik ausgeben ===
% % %     tooLowBins = find(binCounts < minPerBin);
% % %     tooHighBins = find(binCounts > maxPerBin);
% % % 
% % %     if ~isempty(tooLowBins)
% % %         warning('⚠️  %d Bins sind unterbesetzt (< %d Events): %s', ...
% % %                 length(tooLowBins), minPerBin, mat2str(tooLowBins));
% % %     end
% % %     if ~isempty(tooHighBins)
% % %         warning('⚠️  %d Bins sind überbesetzt (> %d Events): %s', ...
% % %                 length(tooHighBins), maxPerBin, mat2str(tooHighBins));
% % %     end
% % % 
% % %     % Optional: Binverteilung als Tabelle anzeigen
% % %     fprintf('\n=== Bin-Auslastung (min %d, ideal %d, max %d) ===\n', ...
% % %         minPerBin, idealPerBin, maxPerBin);
% % %     disp(table((1:nTimeBins)', binCounts', ...
% % %          'VariableNames', {'Bin', 'EventCount'}));
% % % end
% % % 
% % % % 
% % % % % function event_onsets = generate_event_onsets_binned(nTrialsPerCondition, frameRate, nTimeBins, flag_plot)
% % % % % % GENERATE_EVENT_ONSETS_BINNED
% % % % % % Events werden über Zeitbins möglichst gleich verteilt, aber mit 800ms Mindestabstand.
% % % % % %
% % % % % % INPUT:
% % % % % %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–4)
% % % % % %   frameRate           – Hz (z. B. 120)
% % % % % %   nTimeBins           – Anzahl gewünschter Zeitbins
% % % % % %   flag_plot           – true/false für Visualisierung
% % % % % %
% % % % % % OUTPUT:
% % % % % %   event_onsets – Zell-Array {1x4} mit je [nTrials x maxEvents] Eventframe-Matrizen
% % % % % 
% % % % %     rng('default');
% % % % %     rng(1);
% % % % %     maxTries = 10000;
% % % % %     binMargin = 0.05;  % 10% Toleranz bei Zielanzahl pro Bin
% % % % % 
% % % % %     precuetime = 1.7;
% % % % %     totalDuration = 3.0;
% % % % %     totalFrames = round(totalDuration * frameRate);
% % % % %     minSeparationFrames = round(0.8 * frameRate);
% % % % %     possibleFrames = 1:totalFrames;
% % % % %     timeVec = (possibleFrames - 1) / frameRate - precuetime;
% % % % % 
% % % % %     totalEvents = nTrialsPerCondition * sum(1:4);
% % % % %     idealPerBin = round(totalEvents / nTimeBins);
% % % % %     maxPerBin = ceil(idealPerBin * (1 + binMargin));
% % % % % 
% % % % %     % === Zeitbins vorbereiten ===
% % % % %     binEdges = linspace(timeVec(1), timeVec(end), nTimeBins + 1);
% % % % %     binCounts = zeros(1, nTimeBins);
% % % % % 
% % % % %     % Für schnellen Frame → Bin Mapping
% % % % %     binIndexPerFrame = discretize(timeVec, binEdges);
% % % % % 
% % % % %     % === Trialweise generieren ===
% % % % %     event_onsets = cell(1, 4);
% % % % %     usedFrames = false(1, totalFrames);
% % % % %     allEventFrames = [];
% % % % % 
% % % % %     for nEvents = 4:-1:1
% % % % %         trials = zeros(nTrialsPerCondition, 4);
% % % % %         for t = 1:nTrialsPerCondition
% % % % %             success = false;
% % % % %             tries = 0;
% % % % %             while ~success && tries < maxTries
% % % % %                 tries = tries + 1;
% % % % % 
% % % % %                 % Erstelle Framepool: nur aus Bins mit unterdurchschnittlicher Belegung
% % % % %                 eligibleBins = find(binCounts < maxPerBin);
% % % % %                 eligibleFrames = find(ismember(binIndexPerFrame, eligibleBins) & ~usedFrames);
% % % % % 
% % % % %                 if length(eligibleFrames) < nEvents
% % % % %                     warning('Nicht genug Frames verfügbar für %d Events in Trial %d', nEvents, t);
% % % % %                     break;
% % % % %                 end
% % % % % 
% % % % %                 % Wähle Kandidaten
% % % % %                 candidate = sort(randsample(eligibleFrames, nEvents));
% % % % % 
% % % % %                 % Prüfe Mindestabstand
% % % % %                 if all(diff(candidate) >= minSeparationFrames)
% % % % %                     % Akzeptiere Trial
% % % % %                     trials(t, 1:nEvents) = candidate;
% % % % %                     usedFrames(candidate) = true;
% % % % %                     allEventFrames = [allEventFrames, candidate];
% % % % % 
% % % % %                     % Aktualisiere Bin-Zähler
% % % % %                     binsUsed = binIndexPerFrame(candidate);
% % % % %                     for b = binsUsed
% % % % %                         if ~isnan(b)
% % % % %                             binCounts(b) = binCounts(b) + 1;
% % % % %                         end
% % % % %                     end
% % % % % 
% % % % %                     success = true;
% % % % %                 end
% % % % %             end
% % % % % 
% % % % %             if ~success
% % % % %                 warning('Trial %d mit %d Events konnte nach %d Versuchen nicht erstellt werden.', ...
% % % % %                         t, nEvents, maxTries);
% % % % %                 % Trials trotzdem nehmen aber warnen
% % % % %                 trials(t, 1:nEvents) = candidate;
% % % % %                 usedFrames(candidate) = true;
% % % % %                 allEventFrames = [allEventFrames, candidate];
% % % % %             end
% % % % %         end
% % % % %         event_onsets{nEvents} = trials;
% % % % %     end
% % % % % 
% % % % %     % === Visualisierung ===
% % % % %     if flag_plot
% % % % %         figure;
% % % % %         colors = lines(4);
% % % % %         allTimes = [];
% % % % % 
% % % % %         for nEvents = 1:4
% % % % %             subplot(3, 2, nEvents);
% % % % %             onsetFrames = event_onsets{nEvents};
% % % % %             onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
% % % % %             for i = 1:nTrialsPerCondition
% % % % %                 x = onsetTimesSec(i, 1:nEvents);
% % % % %                 y = i * ones(1, nEvents);
% % % % %                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
% % % % %                 allTimes = [allTimes, x];
% % % % %             end
% % % % %             title(sprintf('%d Event(s) pro Trial', nEvents));
% % % % %             xlabel('Zeit (s)');
% % % % %             ylabel('Trial');
% % % % %             xlim([min(timeVec), max(timeVec)]);
% % % % %             ylim([0 nTrialsPerCondition + 1]);
% % % % %             grid on;
% % % % %         end
% % % % % 
% % % % %         % Globale Histogramm-Verteilung
% % % % %         subplot(3, 2, [5 6]);
% % % % %         histogram(allTimes, 'BinEdges', binEdges, ...
% % % % %                   'FaceColor', [0.2 0.4 0.8], 'Normalization', 'count');
% % % % %         hold on;
% % % % %         yline(idealPerBin, '--r', 'Ziel');
% % % % %         yline(maxPerBin, ':r', 'Toleranzgrenze');
% % % % %         xlabel('Zeit (s)');
% % % % %         ylabel('Eventanzahl pro Bin');
% % % % %         title('Globale Verteilung über Zeitbins');
% % % % %         subtitle(sprintf('Bin: %2.0f ms', 1000*min(diff(binEdges))))
% % % % %         xlim([min(timeVec), max(timeVec)]);
% % % % %         grid on;
% % % % % 
% % % % %         sgtitle(sprintf('Events mit Mindestabstand (800 ms) über %d Zeitbins', nTimeBins));
% % % % %     end
% % % % % end
