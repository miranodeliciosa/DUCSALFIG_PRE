function event_onsets = generate_event_onsets(nTrialsPerCondition, frameRate, flag_plot, targetFunction)
% GENERATE_EVENT_ONSETS: Eventonsets nach benutzerdefinierter Zeitverteilung
%
% INPUTS:
%   nTrialsPerCondition – Trials pro Eventanzahl (1–4)
%   frameRate           – Hz (z. B. 120)
%   flag_plot           – true/false für Visualisierung
%   targetFunction      – @(t) ... Funktion, die gewünschte Wahrscheinlichkeitsdichte angibt
%
% OUTPUT:
%   event_onsets – Zell-Array mit je [nTrials x maxEvents] Eventframe-Matrizen

    totalDuration = 3.5;
    totalFrames = round(totalDuration * frameRate);
    minSeparationFrames = round(0.8 * frameRate);
    possibleFrames = 1:totalFrames;
    precuetime = 1.7;
    timeVec = (possibleFrames - 1) / frameRate - precuetime;
    numEvents = 3;

    rng('default');
    rng(3); % Reproduzierbar
    event_onsets = cell(1, numEvents);
    allEventFrames = [];

    % === Zielverteilung als Wahrscheinlichkeiten berechnen ===
    targetWeights = targetFunction(timeVec);
    targetWeights = max(targetWeights, 0);           % Negativwerte verhindern
    targetWeights = targetWeights / sum(targetWeights); % Normalisieren auf 1

    maxTries = 10000; % Absicherung

    for nEvents = numEvents:-1:1
        trials = zeros(nTrialsPerCondition, numEvents);
        for t = 1:nTrialsPerCondition
            success = false;
            tries = 0;
            while ~success && tries < maxTries
                tries = tries + 1;
                candidate = sort(randsample(possibleFrames, nEvents, true, targetWeights));
                if all(diff(candidate) >= minSeparationFrames)
                    trials(t, 1:nEvents) = candidate;
                    allEventFrames = [allEventFrames, candidate];
                    success = true;
                end
            end
            if ~success
                warning('Trial %d mit %d Events konnte nicht erzeugt werden.', t, nEvents);         
            end
        end
        event_onsets{nEvents} = trials;
    end

    % === Plot (optional)
    if flag_plot
        figure;
        colors = lines(4);

        for nEvents = 1:numEvents
            subplot(3, 2, nEvents);
            onsetFrames = event_onsets{nEvents};
            onsetTimesSec = (onsetFrames - 1) / frameRate - precuetime;
            for i = 1:nTrialsPerCondition
                x = onsetTimesSec(i, 1:nEvents);
                y = i * ones(1, nEvents);
                scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
            end
            title(sprintf('%d Event(s) pro Trial', nEvents));
            xlabel('Zeit (s)');
            ylabel('Trial');
            xlim([-precuetime 1.5]);
            ylim([0 nTrialsPerCondition + 1]);
            grid on;
        end

        % Histogramm der tatsächlichen Verteilung + Zielverteilung
        subplot(3,2,[5 6]);
        eventTimes = (allEventFrames - 1) / frameRate - precuetime;
       % histogram(eventTimes, 'BinEdges', timeVec, 'Normalization', 'pdf', ...
       %           'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'none'); hold on;
        h = histogram(eventTimes, 'BinWidth', 30/frameRate,...
                  'FaceColor', [0.2 0.2 0.8]); hold on;%
        plot(timeVec, targetWeights * frameRate * max(h.Values), 'r-', 'LineWidth', 2);
        xlabel('Zeit (s)');
        ylabel('Event-Anzahl');
        legend('Resultierende Verteilung', 'Zielverteilung');
        title('Globale Eventverteilung (empirisch vs. Ziel)');
        xlim([-precuetime 1.5]);
        grid on;

        sgtitle('Event-Onsets nach Zielzeitverteilung');
    end
end



% %VERSION 4
% function event_onsets = generate_event_onsets(nTrialsPerCondition, frameRate, flag_plot)
% % GENERATE_EVENT_ONSETS: gleichverteilte, abstandsgeprüfte Eventonsets
% %
% % INPUTS:
% %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–4)
% %   frameRate           – Hz (z. B. 120)
% %   flag_plot           – true/false für Visualisierung
% %
% % OUTPUT:
% %   event_onsets – Zell-Array mit je [nTrials x maxEvents] Eventframe-Matrizen
% 
%     totalDuration = 3.5;
%     totalFrames = round(totalDuration * frameRate);
%     minSeparationFrames = round(0.8 * frameRate);
%     possibleFrames = 1:totalFrames;
%     preCueTime = 2;
%     timeVec = (possibleFrames - 1) / frameRate - preCueTime;
% 
%     rng('default')
%     rng(3); % reproduzierbar
%     event_onsets = cell(1, 4);
%     allEventFrames = [];
%     frameUsage = zeros(1, totalFrames);  % globaler Frame-Zähler
% 
%     maxTries = 10000;  % absichern gegen Endlosschleifen
% 
%     for nEvents = 1:4
%         trials = zeros(nTrialsPerCondition, 4);  % max 4 Events pro Trial
%         for t = 1:nTrialsPerCondition
%             success = false;
%             tries = 0;
%             while ~success && tries < maxTries
%                 tries = tries + 1;
% 
%                 % === Gewichtete Samplingwahrscheinlichkeit (Inverse der Nutzung)
%                 weights = max(max(frameUsage) - frameUsage + 1, 1);
%                 weights = weights / sum(weights);
% 
%                 % Ziehe Kandidaten
%                 candidate = sort(randsample(possibleFrames, nEvents, true, weights));
% 
%                 % Prüfe Mindestabstand zwischen Events im Trial
%                 if all(diff(candidate) >= minSeparationFrames)
%                     trials(t, 1:nEvents) = candidate;
%                     frameUsage(candidate) = frameUsage(candidate) + 1;
%                     allEventFrames = [allEventFrames, candidate];
%                     success = true;
%                 end
%             end
% 
%             if ~success
%                 warning('Trial %d mit %d Events konnte nach %d Versuchen nicht erzeugt werden.', t, nEvents, maxTries);
%             end
%         end
%         event_onsets{nEvents} = trials;
%     end
% 
%     % === Plot (optional)
%     if flag_plot
%         figure;
%         colors = lines(4);
%         for nEvents = 1:4
%             subplot(3, 2, nEvents);
%             onsetFrames = event_onsets{nEvents};
%             onsetTimesSec = (onsetFrames - 1) / frameRate - preCueTime;
%             for i = 1:nTrialsPerCondition
%                 x = onsetTimesSec(i, 1:nEvents);
%                 y = i * ones(1, nEvents);
%                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
%             end
%             title(sprintf('%d Event(s) pro Trial', nEvents));
%             xlabel('Zeit (s)');
%             ylabel('Trial');
%             xlim([min(timeVec) max(timeVec)]);
%             ylim([0 nTrialsPerCondition + 1]);
%             grid on;
%         end
% 
%         % Histogramm der globalen Verteilung
%         subplot(3,2,[5 6]);
%         eventTimes = (allEventFrames - 1) / frameRate - preCueTime;
%         histogram(eventTimes, 'BinWidth', 30/frameRate, 'FaceColor', [0.2 0.2 0.8]);
%         xlabel('Zeit (s)');
%         ylabel('Event-Anzahl');
%         title('Globale Verteilung aller Event-Onsets');
%         xlim([min(timeVec) max(timeVec)]);
%         grid on;
% 
%         sgtitle('Gleichverteilte Event-Onsets mit Mindestabstand');
%     end
% end



% VERSION 2
% function event_onsets = generate_event_onsets(nTrialsPerCondition, frameRate, flag_plot, method)
% % GENERATE_EVENT_ONSETS erstellt Onsets für 1–4 Events pro Trial
% % mit Mindestabstand von 800 ms, optionaler Visualisierung und kontrollierter Verteilung.
% %
% % INPUTS:
% %   nTrialsPerCondition – Anzahl der Trials pro Eventanzahl (1–4)
% %   frameRate           – Bildschirmwiederholrate in Hz (z. B. 120)
% %   flag_plot           – true/false, ob geplottet werden soll
% %   method              – 'random' (default) oder 'balanced'
% %
% % OUTPUT:
% %   event_onsets – Zell-Array {1x4} mit je [nTrials x maxEvents] Eventframe-Matrizen
% 
%     if nargin < 4
%         method = 'random';
%     end
% 
%     totalDuration = 3.0; % Sekunden (-1.5 bis +1.5)
%     totalFrames = round(totalDuration * frameRate);
%     minSeparationFrames = round(0.8 * frameRate);
%     possibleFrames = 1:totalFrames;
%     event_onsets = cell(1, 4);
%     allEventFrames = [];
%     frameUsage = zeros(1, totalFrames);  % für 'balanced'-Modus
% 
%     rng('default');
%     rng(1); % Reproduzierbarkeit
% 
%     for nEvents = 1:4
%         trials = zeros(nTrialsPerCondition, 4);
%         for t = 1:nTrialsPerCondition
%             success = false;
%             while ~success
%                 switch method
%                     case 'random'
%                         candidate = sort(randsample(possibleFrames, nEvents));
%                     case 'balanced'
%                         % Gewichtete Auswahl: bevorzuge selten benutzte Frames
%                         weights = max(max(frameUsage) - frameUsage + 1, 1);
%                         weights = weights / sum(weights);
%                         candidate = sort(randsample(possibleFrames, nEvents, true, weights));
%                     otherwise
%                         error('Ungültige Methode. Verwende "random" oder "balanced".');
%                 end
% 
%                 % Mindestabstand prüfen
%                 if all(diff(candidate) >= minSeparationFrames)
%                     trials(t, 1:nEvents) = candidate;
%                     success = true;
%                     frameUsage(candidate) = frameUsage(candidate) + 1;
%                     allEventFrames = [allEventFrames, candidate];
%                 end
%             end
%         end
%         event_onsets{nEvents} = trials;
%     end
% 
%     % ===== VISUALISIERUNG =====
%     if flag_plot
%         figure;
%         colors = lines(4);
% 
%         % Subplots für 1–4 Events
%         for nEvents = 1:4
%             subplot(3, 2, nEvents);
%             onsetFrames = event_onsets{nEvents};
%             onsetTimesSec = (onsetFrames - 1) / frameRate - 1.5;
%             for i = 1:nTrialsPerCondition
%                 x = onsetTimesSec(i, 1:nEvents);
%                 y = i * ones(1, nEvents);
%                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
%             end
%             title(sprintf('%d Event(s) pro Trial', nEvents));
%             xlabel('Zeit (s)');
%             ylabel('Trial');
%             xlim([-1.5 1.5]);
%             ylim([0 nTrialsPerCondition+1]);
%             grid on;
%         end
% 
%         % Globales Histogramm aller Events
%         subplot(3,2,[5 6]);
%         eventTimes = (allEventFrames - 1) / frameRate - 1.5;
%         histogram(eventTimes, 'BinWidth', 1/frameRate, 'FaceColor', [0.2 0.2 0.8]);
%         xlabel('Zeit (s)');
%         ylabel('Event-Anzahl');
%         title(sprintf('Gesamte Eventverteilung (%s)', method));
%         xlim([-1.5 1.5]);
%         grid on;
% 
%         sgtitle('Event-Onsets und zeitliche Verteilung');
%     end
% end


% VERSION 1
% function event_onsets = generate_event_onsets(nTrialsPerCondition, frameRate, flag_plot)
% % GENERATE_EVENT_ONSETS erstellt Onsets für 1–4 Events pro Trial
% % unter Einhaltung eines Mindestabstands von 800 ms und optionaler Visualisierung.
% %
% % INPUTS:
% %   nTrialsPerCondition – Anzahl der Trials für jede Eventanzahl (1–4)
% %   frameRate – Bildschirmwiederholrate in Hz (z. B. 120)
% %   flag_plot – true/false, ob die Visualisierung angezeigt werden soll
% %
% % OUTPUT:
% %   event_onsets – Zelle mit 4 Zellen (1–4 Events), jede enthält eine 
% %                  [nTrialsPerCondition x maxEvents] Matrix mit Frame-Onsets
% 
%     totalDuration = 3.5; % Sekunden (-1.5s bis +1.5s)
%     totalFrames = round(totalDuration * frameRate); % z. B. 360 bei 120Hz
%     minSeparationFrames = round(0.8 * frameRate);   % z. B. 96 Frames
% 
%     timeVec = linspace(-2, 1.5, totalFrames);
%     possibleFrames = 1:totalFrames;
% 
%     rng('default');
%     rng(2); % für Reproduzierbarkeit
%     event_onsets = cell(1, 4); % für 1–4 Events
% 
%     % Sammle alle Event-Frames für Histogramm
%     allEventFrames = [];
% 
%     for nEvents = 1:4
%         trials = zeros(nTrialsPerCondition, 4); % Platz für max. 4 Events
%         for t = 1:nTrialsPerCondition
%             success = false;
%             while ~success
%                 candidate = sort(randsample(possibleFrames, nEvents));
%                 if all(diff(candidate) >= minSeparationFrames)
%                     trials(t, 1:nEvents) = candidate;
%                     success = true;
%                 end
%             end
%         end
%         event_onsets{nEvents} = trials;
%         allEventFrames = [allEventFrames trials(:, 1:nEvents)];
%     end
% 
%     % === Optionale Visualisierung ===
%     if flag_plot
%         figure;
%         colors = lines(4);
% 
%         % Plot 1–4 Event-Trials
%         for nEvents = 1:4
%             subplot(3, 2, nEvents);
%             onsetFrames = event_onsets{nEvents};
%             % Konvertiere in Sekunden für Plot
%             onsetTimesSec = (onsetFrames - 1) / frameRate - 1.5;
%             for i = 1:nTrialsPerCondition
%                 y = i * ones(1, nEvents);
%                 x = onsetTimesSec(i, 1:nEvents);
%                 scatter(x, y, 20, colors(nEvents,:), 'filled'); hold on;
%             end
%             title(sprintf('%d Event(s) pro Trial', nEvents));
%             xlabel('Zeit (s)');
%             ylabel('Trial');
%             xlim([-1.5 1.5]);
%             ylim([0 nTrialsPerCondition+1]);
%             grid on;
%         end
% 
%         % Plot globales Histogramm der Eventhäufigkeit über Zeit
%         subplot(3,2,[5 6]);
%         % Konvertiere Event-Frames zu Zeit
%         eventTimes = (allEventFrames(:) - 1) / frameRate - 1.5;
%         histogram(eventTimes, 'BinWidth', 10/frameRate, 'FaceColor', [0.2 0.2 0.8]);
%         xlabel('Zeit (s)');
%         ylabel('Anzahl Events');
%         title('Gesamte Eventverteilung über Zeit');
%         xlim([-1.5 1.5]);
%         grid on;
% 
%         sgtitle('Event-Onsets und globale Verteilung');
%     end
% end
