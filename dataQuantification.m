% To be done after PostCuration_BK.m
% General data input: single units (array) stimuli order (table), aligned spikes (cell array),
% Post processing data steps:
%   1. get firing rate during set time frame for each trial of each unit
%   2. quatify which cluster responds to which event
%       responsive unit = sig different firing after stim vs baseline
%   3. quantify response strength
%   4. combine data from multiple animals

clearvars

%SxA, SOM session matching 
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % M9

%% ----------------------- FRA analysis & Plotting ----------------------- %%
% output: FRA & MedFSL 4D: intensity, frequency, set number, cluster

close all

OutPath = 'D:\DATA\Processed\M26\postCA'; % output directory
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M26\'; % stimuli parameters

% load FRA session(s) aligned spikes
aligned_spikes_files = dir(fullfile(OutPath, '*FRA_AlignedSpikes.mat'));

for file = 1:size(aligned_spikes_files, 1)

    % load each FRA aligned spikes file
    %aligned_spikes = load([aligned_spikes_files(file).folder '\' aligned_spikes_files(file).name]);
    % load corresponding stimuli file
    %session = aligned_spikes_files(file).name(6:7);
    %stim_files = dir(fullfile(BehaviorPath, ['*_S' session '_*.mat']));
    %stimuli_parameters = load([stim_files.folder '\' stim_files.name]);

    session = str2double(aligned_spikes_files(file).name(6:7));
    [cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);

    % FRA analysis saves heatmap figures
    FSL = 0;
    heatmap = 1;
    FRAanalysis(stimuli_parameters, aligned_spikes, cids, OutPath, heatmap, FSL);

end

%% quantitative analysis
clearvars
%set paths
mousenr = 24;
OutPath = 'D:\DATA\Processed\M24';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M24';

%sessions = [1 2 3 4 5 6 7 8 9 10];% M10
%sessions = [1 2 3 4 5 6 7 8 9 10];% M11
%sessions = [1 2 6 7 8 9 10 11 12 13 14 15 16]; % M19 all sessions to be analyzed 
%sessions = [1 2 3 4 5 6 7 8 9 10 11 12];% M20
sessions = [15, 16];

for session = 1:2 %1:length(sessions)
    if isempty(dir([OutPath '\*_S' num2str(sessions(session),'%02i') '*_ResponseProperties.mat']))
        dataQuantification_analysis(mousenr, OutPath, BehaviorPath, sessions(session))
    else
        continue
    end
end
%% Pool datasets across mice
stimOrder = readtable('D:\DATA\Processed\M10-11-19-20_stimOrder.csv');
stimOrder.AM_pre(2) = NaN; % discard M11 AM session, wrong params
fn = 'M10-11-19-20';

dataQuantification_poolDatasets(stimOrder, fn)


%% TO DO: cycle analysis
% define analysis window - FR/vibrotactile cycle
% analysis window dependent on cylce
% number of data points dependent on cycle

OutPath = 'D:\DATA\Processed\M20\ICX';
BehaviorPath = 'D:\DATA\Behavioral Stimuli\M20';

% load data
[~, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, 3, BehaviorPath);
freqs = unique(stimuli_parameters.Stm.SomFreq);
% max size matrix
minwindow = max(freqs)/1000;
minwindow * max(stimuli_parameters.Stm.SomDur)/1000;

if strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine') % && strcmp(stimuli_parameters.Par.Rec, "SxA")
    for freq = 1:length(freqs)
        PreT = (str2double(stimuli_parameters.Par.SomatosensoryISI)/3)/1000; % baseline period
        PostT = str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;

        firingrate(aligned_spikes, PreT, PostT)
    end
end

%% CYCLE QUANTIFICATION - IN PROGRESS

OutPath = 'V:\Data\ProcessedData\Blom\Processed\M20\ICX';
BehaviorPath = 'V:\Data\InVivoEphys\Blom\BehavioralStimuli\M20';
% OutPath = 'V:\Data\ProcessedData\Blom\Processed\M19\ICX';
% BehaviorPath = 'V:\Data\InVivoEphys\Blom\BehavioralStimuli\M19';

% load data
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, ~] = loadData(OutPath, 3, BehaviorPath);
%%
cluster = 11;
spikeTimes = aligned_spikes(:,cluster);

SomFreq = 50;
Amplitude = 0.3;
Delay=18e-3;
for MMType = ["SO","SA"]
SomDelay = 0.25;
s_idx =     stimuli_parameters.Stm.SomFreq == SomFreq & ...
            stimuli_parameters.Stm.Amplitude == Amplitude & ...
            strcmp(stimuli_parameters.Stm.MMType, MMType);

SpkT = spikeTimes(s_idx);
if strcmp(MMType,'SA')
    for ii = 1:length(SpkT)
        SpkT{ii} = SpkT{ii} - SomDelay;
    end
end
stim = stimuli_parameters.Stm(s_idx,:);
stim = stim(1,:);
StimDur = 0.001*stim.SomDur;

    % prepare for display
    switch MMType
        case "SA"
            fig = figure(2);
        case "SO"
            fig = figure(1);
        otherwise
            fig = figure(3);
    end
    clf(fig)
    xRange = [-0.3,StimDur+0.3];

% 1) get cycle times

Mf = SomFreq;
[CycT,edges] = CycTimesPerCycle(SpkT,StimDur, Mf,Delay);

% 2) Sum the number of elements in each column of CycT
NStim = size(CycT,1);
NCyc = size(CycT,2);
NumSpikesPerCyc = cellfun(@numel, CycT);
avgNumSpikesPerCyc = mean(NumSpikesPerCyc,1);
    % display
    ncols = 5;
    % plot stimlus
    ax0 = subplot(ncols,1,1);
    tt = 0:0.0001:StimDur;
    stim_waveform = (1-cos(2*pi*SomFreq.*tt));
    plot(ax0,tt,stim_waveform,'k-')
    xlim(ax0,xRange)
    

    ax1 = subplot(ncols,1,[2 3]);
    binCenters = 0.5.*(edges(1:end-1) + edges(2:end));
    hold off
    
    % plot bin count
    yyaxis(ax1,'left')
    plot(ax1,binCenters,NStim.*avgNumSpikesPerCyc,'x-'); hold on;
    xline(ax1,edges,'--','Color',[.7,.7,.7])
    xlim(ax1,xRange)
    period = 1/SomFreq;

    % plot PSTH
    hist_edges = xRange(1):(min(0.1*period,0.005)):(max(StimDur,xRange(2)));
    histogram(ax1,vertcat(SpkT{:}),hist_edges)
    ylabel(ax1,'Number of spikes')


% 3) latency
meanLatencyPerCyc   =   nan(1,NCyc);
stdLatencyPerCyc    =   nan(1,NCyc);
for cyc_idx = 1:NCyc
    tSpk = vertcat(CycT{:,cyc_idx});
    meanLatencyPerCyc(cyc_idx) = mean(tSpk);
    stdLatencyPerCyc(cyc_idx) = std(tSpk);
end
    % display
    hold(ax1,'on');
    yyaxis(ax1,'right')
    plot(ax1,meanLatencyPerCyc,stdLatencyPerCyc*1000,'s:k')
    ylabel(ax1,'Latency jitter (ms)')
    ylim(ax1,[0,inf])
    hold(ax1,'off');
   
% 4) Vector strength
% Normalize all values in CycT: subtract edges(m) and divide by period
period = 1 / SomFreq;
CycT_norm = CycT;
for n = 1:NStim
    for m = 1:NCyc
        CycT_norm{n,m} = (CycT{n,m} - edges(m)) / period;
    end
end
VSPerCycle = nan(1,NCyc);
PhasePerCycle = nan(1,NCyc);
RayleighPerCycle = nan(1,NCyc);
for m = 1:NCyc
    [VSPerCycle(m),PhasePerCycle(m),RayleighPerCycle(m)] = calcVS(CycT_norm(:,m));
end
    % display
    alpha = 0.05;
    ax2 = subplot(ncols,1,[4 5]);
    yyaxis(ax2,'left')
    hold(ax2,"off")
    plot(ax2,binCenters,VSPerCycle,'o-')
    hold(ax2,"on")
    scatter(ax2,binCenters(RayleighPerCycle<alpha),VSPerCycle(RayleighPerCycle<alpha),'o','filled')
    hold(ax2,"off")
    ylabel("Vector strength")
    ylim(ax2,[0,1])

    yyaxis(ax2,'right')
    plot(ax2,binCenters,PhasePerCycle)
    ylabel("Delay-subtracted phase (cycle)")
    ylim(ax2,[0,1])
    xline(ax2,edges,'--','Color',[.7,.7,.7])
    xlim(ax2,xRange)

    sgtitle(['Cluster ', num2str(cids(cluster)), ...
            newline,...
            '    Condition: ',char(MMType), ...
            '    Freq: ', num2str(SomFreq,'%d Hz'), ...
            '    Amp: ', num2str(Amplitude,'%.1f V'), ...
            newline,...
            '    Analysis delay: ',num2str(Delay*1000,'%d ms')])
    linkaxes([ax0,ax1,ax2],'x')
end
%% Checking group delay via phase
cluster = 11;
spikeTimes = aligned_spikes(:,cluster);

Amplitude = 0.3;
SomFreqList = [10,20,50,100];
MMType = 'SO';
Delay = 0;
DelayList = 0:1e-3:30e-3;
meanPhaseCorr = nan(length(DelayList),length(SomFreqList));
for d_idx = 1:length(DelayList)
    Delay = DelayList(d_idx);
    meanPhase = nan(size(SomFreqList));
    for f_idx = 1:length(SomFreqList)
        SomFreq = SomFreqList(f_idx);
        s_idx =     stimuli_parameters.Stm.SomFreq == SomFreq & ...
                    stimuli_parameters.Stm.Amplitude == Amplitude & ...
                    strcmp(stimuli_parameters.Stm.MMType, MMType);
        
        SpkT = spikeTimes(s_idx);
        
        [CycT,edges] = CycTimesPerCycle(SpkT,StimDur,SomFreq,Delay);
        NStim = size(CycT,1);
        NCyc = size(CycT,2);
    
        period = 1 / SomFreq;
        CycT_norm = CycT;
        for n = 1:NStim
            for m = 1:NCyc
                CycT_norm{n,m} = (CycT{n,m} - edges(m)) / period;
            end
        end
        VSPerCycle = nan(1,NCyc);
        PhasePerCycle = nan(1,NCyc);
        RayleighPerCycle = nan(1,NCyc);
        for m = 1:NCyc
            [VSPerCycle(m),PhasePerCycle(m),RayleighPerCycle(m)] = calcVS(CycT_norm(:,m));
        end
        meanPhase(f_idx) = mean(PhasePerCycle(1:3),'omitnan');
    end
    meanPhaseCorr(d_idx,:) = meanPhase;
end
figure(3);
% meanPhaseCorr = mod(meanPhase - DelayList' ./ (1./SomFreqList),1);
plot(DelayList,meanPhaseCorr)
legend(string(SomFreqList))
ylim([0,1])
%%
% Example cell array structure: rows (trials), columns (neurons)
spikeTimes = aligned_spikes;

% Define fixed stimulus duration (in seconds)
stimulusDuration = max(stimuli_parameters.Stm.SomDur)/1000; % 500 ms

% Initialize result cell array
[nTrials, clusters] = size(spikeTimes);
spikeCountsPerCycle = cell(nTrials, clusters);

% Loop through neurons first
for cluster = 1:clusters
    for trial = 1:nTrials
        freq = stimuli_parameters.Stm.SomFreq(trial); % Get stimulus frequency per trial

        % Define cycle edges based on stimulus duration
        if freq == 0 % one bin for full stimulus
            cycleEdges = [0, stimulusDuration]; 
        else
            cycleDuration = 1 / freq; % cycle duration
            cycleEdges = 0:cycleDuration:stimulusDuration;
        end

        spikes = spikeTimes{trial, cluster}; % extract spike times
        [N, edges] = histcounts(vertcat(aligned_spikes{SA, cluster}), preT:binsize:postT); % SA

        if ~isempty(spikes) % Handle empty spike cases
            spikeCountsPerCycle{trial, cluster} = histcounts(spikes, cycleEdges); % Count spikes per cycle
            % histcounts(vertcat(aligned_spikes{SA, cluster}), preT:binsize:postT); % SA
        else
            spikeCountsPerCycle{trial, cluster} = zeros(1, length(cycleEdges) - 1); % Fill with zeros if no spikes
        end
    end
end

close all

ufreqs = unique(stimuli_parameters.Stm.SomFreq);

for cluster = 1:clusters
    % figure per cluster with all conditions
    figure;
    numSubplots = length(ufreqs);

    for freqIdx = 1:length(ufreqs)
        freq = ufreqs(freqIdx);
        subplot(numSubplots, 1, freqIdx); % Arrange subplots in a vertical layout
        hold on;

        % Find trials for this frequency
        %trialIndices = find(stimuli_parameters.Stm.SomFreq == freq & stimuli_parameters.Stm.Var25 == 3);
        if freq == 0
            %continue
            trialIndices = find((stimuli_parameters.Stm.Var25 == 1));
        else
            trialIndices = find((stimuli_parameters.Stm.Amplitude ~= 0.1) & (stimuli_parameters.Stm.SomFreq == freq) & (stimuli_parameters.Stm.Var25 == 3));
        end

        % Loop through trials of the current stimulus frequency
        for i = 1:length(trialIndices)
            trialIdx = trialIndices(i);
            spikeCounts(:,1) = spikeCountsPerCycle{trialIdx, cluster};
            plot(1:length(spikeCounts), spikeCounts, '-o');

        end

        %plot(1:length(spikeCounts), sum(spikeCountsPerCycle{trialIndices, cluster}), 'k-o');

        xlabel('Cycle Index');
        ylabel('Spike Count');
        title(sprintf('%d Hz', freq));

    end

    sgtitle(['Cluster ' num2str(cids(cluster))])
    hold off;

end



%% Example cell array containing spike times (in seconds)
spikeTimes = { [0.02, 0.04, 0.07, 0.12], [0.03, 0.06, 0.09, 0.15] }; 

% Define stimulus frequency (in Hz)
freq = 50; 
cycleDuration = 1 / freq;  % Duration of one cycle in seconds

% Initialize result array
spikeCountsPerCycle = cell(size(spikeTimes));

% Loop through each cell in the array
for i = 1:numel(spikeTimes)
    spikes = spikeTimes{i}; % Extract spike times for this trial
    cycleEdges = 0:cycleDuration:max(spikes) + cycleDuration; % Define cycle boundaries
    spikeCountsPerCycle{i} = histcounts(spikes, cycleEdges); % Count spikes per cycle
end

% Display results
disp(spikeCountsPerCycle);

%% STATS - to do

%% OPTIONAL: add SOM responses to SxA table, if needed
mousenr = 9; % M11
%session = [3, 4]; % M11
%session = [2, 5]; % M10
%session = [5, 3]; % M8
%session = [4, 6]; % SxA, SOM M9

clearvars

OutPath = 'D:\DATA\Processed\M10-11-19-20';

% load SxA data
load('D:\DATA\Processed\M10-11-19-20\M10-11-19-20_PxA_pre_UnitResponses.mat')
StimResponseFiring_all_SxA = StimResponseFiring_all;
unitResponses_SxA = StimResponseFiring_all.unitResponses;
firingmean_SxA = firing_mean; % size: 5     9     4    24
MouseNum_SxA = StimResponseFiring_all.MouseNum;
session_SxA = StimResponseFiring_all.session;
amplitudes_SxA = StimResponseFiring_all.amplitudes;
frequencies_SxA = StimResponseFiring_all.frequencies;
PostT_SxA = StimResponseFiring_all.PostT;
PreT_SxA = StimResponseFiring_all.PreT;
FSLmed_SxA = StimResponseFiring_all.FSLmed;
FSLiqr_SxA = StimResponseFiring_all.FSLiqr;
hvalue_SxA = StimResponseFiring_all.hvalue;

% load pressure data
load('D:\DATA\Processed\M10-11-19-20\M10-11-19-20_SOM_P_pre_UnitResponses.mat')
StimResponseFiring_all_SOM = StimResponseFiring_all;
unitResponses_SOM = StimResponseFiring_all.unitResponses;
firingmean_SOM = firing_mean; % size: 1     9     1    40
MouseNum_SOM = StimResponseFiring_all.MouseNum;
session_SOM = StimResponseFiring_all.session;
amplitudes_SOM = StimResponseFiring_all.amplitudes;
frequencies_SOM = StimResponseFiring_all.frequencies;
PostT_SOM = StimResponseFiring_all.PostT;
PreT_SOM = StimResponseFiring_all.PreT;
FSLmed_SOM = StimResponseFiring_all.FSLmed;
FSLiqr_SOM = StimResponseFiring_all.FSLiqr;
hvalue_SOM = StimResponseFiring_all.hvalue;

% combine data
unitResponses = joindata(unitResponses_SxA, unitResponses_SOM); % combine unit response table

MouseNum = [MouseNum_SOM, MouseNum_SxA]; % combine mousenr
session = [session_SOM, session_SxA]; %combine session nr
amplitudes = [amplitudes_SOM, amplitudes_SxA]; % combine amplitudes

% combine freqs
frequencies_SOM = [-Inf -Inf];
frequencies_SOM(end+1:size(frequencies_SxA,1),:) = missing;
frequencies_SOM(end+1:size(frequencies_SxA,1),:) = missing;
frequencies = [frequencies_SOM, frequencies_SxA]; 

% combine postT and preT
PostT_SOM(end+1:size(PostT_SxA,1),:) = missing;
PreT_SOM(end+1:size(PreT_SxA,1),:) = missing;
PostT = [PostT_SOM, PostT_SxA];
PreT = [PreT_SOM, PreT_SxA];

% combine firing mean data
temp = NaN(size(firingmean_SxA,1),size(firingmean_SOM,2),size(firingmean_SxA,3),size(firingmean_SOM,4)); % make temporary matrix in correct dimensions
temp(1,1,1,:) = firingmean_SOM(1,1,1,:); % add control
temp(1,2:end,3,:) = firingmean_SOM(1,2:end,1,:); % add som only
firingmean_SOM = temp;
firing_mean = cat(4,firingmean_SOM, firingmean_SxA);

% combine FSLmed
temp = NaN(size(FSLmed_SxA,1),size(FSLmed_SOM,2),size(FSLmed_SxA,3),size(FSLmed_SOM,4)); % make temporary matrix in correct dimensions
temp(1,1,1,:) = FSLmed_SOM(1,1,1,:); % add control
temp(1,2:end,3,:) = FSLmed_SOM(1,2:end,1,:); % add som only
FSLmed_SOM = temp;
FSLmed = cat(4,FSLmed_SOM, FSLmed_SxA);

% combine FSLiqr
temp = NaN(size(FSLiqr_SxA,1),size(FSLiqr_SOM,2),size(FSLiqr_SxA,3),size(FSLiqr_SOM,4)); % make temporary matrix in correct dimensions
temp(1,1,1,:) = FSLiqr_SOM(1,1,1,:); % add control
temp(1,2:end,3,:) = FSLiqr_SOM(1,2:end,1,:); % add som only
FSLiqr_SOM = temp;
FSLiqr = cat(4,FSLiqr_SOM, FSLiqr_SxA);

% combine hval
temp = NaN(size(hvalue_SxA,1),size(hvalue_SOM,2),size(hvalue_SxA,3),size(hvalue_SOM,4)); % make temporary matrix in correct dimensions
temp(1,1,1,:) = hvalue_SOM(1,1,1,:); % add control
temp(1,2:end,3,:) = hvalue_SOM(1,2:end,1,:); % add som only
hvalue_SOM = temp;
hvalue = cat(4,hvalue_SOM, hvalue_SxA);

% save in struct
StimResponseFiring_all.StimType = 'PxA_pre';
StimResponseFiring_all.MouseNum = MouseNum;
StimResponseFiring_all.session = session;
StimResponseFiring_all.amplitudes = amplitudes;
StimResponseFiring_all.frequencies = frequencies;
StimResponseFiring_all.PreT = PreT;
StimResponseFiring_all.PostT = PostT;
StimResponseFiring_all.unitResponses = unitResponses;
StimResponseFiring_all.firing_mean = firing_mean;
StimResponseFiring_all.FSLmed = FSLmed;
StimResponseFiring_all.FSLiqr = FSLiqr;
StimResponseFiring_all.hvalue = hvalue;


% include condition when saving
filename = sprintf('M10-11-19-20_PxA_P_pre_UnitResponses');
save(fullfile(OutPath, filename), "firing_mean", "StimResponseFiring_all")

%% OPTIONAL: Bandwith FRA analysis

uAmp = StimResponseFiring_all.amplitudes(:,1); % unique(stimuli_parameters.Stm.Intensity);
uFreq = StimResponseFiring_all.frequencies(:,1); %unique(stimuli_parameters.Stm.Freq);
NClu = length(StimResponseFiring_all.unitResponses.Cluster);  %length(cids);
%NClu = length(selected_cids);
StimResponseFiring = StimResponseFiring_all;

deltaInt = [0,10,20,30,40];

% initiate variables
excBand = NaN(length(deltaInt), 2, NClu);
inhBand = NaN(length(deltaInt), 2, NClu);
excBW_oct = NaN(length(deltaInt),NClu);
inhBW_oct = NaN(length(deltaInt),NClu);
excQ = NaN(length(deltaInt),NClu);
inhQ = NaN(length(deltaInt),NClu);
CF = NaN(NClu,1);
BF = NaN(NClu,1);
BestFreqAmp = NaN(NClu,1);
%
for cluster = 1:NClu

    % only for selected clusters
    if selected_cids_idx(cluster)
    spontRate = StimResponseFiring.firing_mean(1, 1, 1, cluster);

    % find best frequency (BF)
    [~,maxIdx] = max(StimResponseFiring.firing_mean(:, :, 1, cluster),[],'all','linear'); % find max firing rate change
    [maxRow, maxCol] = ind2sub(size(StimResponseFiring.firing_mean(:, :, 1, cluster)), maxIdx); % translate to position in matrix
    if StimResponseFiring.hvalue(maxRow, maxCol, cluster) % select from sig responses
        BF(cluster) = uFreq(maxRow);
        BestFreqAmp(cluster) = uAmp(maxCol);
    end

    % find minimum threshold from significant neurons
    [row, col] = find(StimResponseFiring.hvalue(:,:,cluster)==1); % select from sig responses
    if isempty(col)
        minThr = NaN;
    else
        minThr = uAmp(min(col)); % minimal threshold
    end

    % find minimum threshold for all units
    %[row, col] = find(FRA.FRASR(:,:,1,cluster))

    %find charactristic frequency (CF)
    CF(cluster) = mean(uFreq(row(col == min(col))));  

    % bandwith
    for dInt = 1:length(deltaInt)
        idx = find(uAmp == minThr + deltaInt(dInt)); % min threshold + delta

        excBandIdx = (StimResponseFiring.firing_mean(idx,:,1,cluster) >= spontRate) & (StimResponseFiring.hvalue(idx,:,cluster)==1);
        inhBandIdx = (StimResponseFiring.firing_mean(idx,:,1,cluster) <= spontRate) & (StimResponseFiring.hvalue(idx,:,cluster)==1);

        if(sum(excBandIdx)>0) % excitatory band
            band_start = find(excBandIdx,1,'first');
            excBand(dInt,1, cluster) = sqrt(uFreq(band_start) * uFreq(band_start - 1)) ;
            band_end = find(excBandIdx,1,'last');
            excBand(dInt,2, cluster) = sqrt(uFreq(band_end) * uFreq(band_end + 1)) ;

            excBW_oct(dInt, cluster) = log2(excBand(dInt,2,cluster) / excBand(dInt,1,cluster));
            excQ(dInt, cluster) = CF(cluster) / (excBand(dInt,2, cluster) - excBand(dInt,1, cluster));
        end

        if(sum(inhBandIdx)>0) % inhibitory band
            band_start = find(inhBandIdx,1,'first');
            inhBand(dInt, 1, cluster) = sqrt(uFreq(band_start) * uFreq(band_start - 1)) ;
            band_end = find(inhBandIdx, 1,'last');
            inhBand(dInt, 2, cluster) = sqrt(uFreq(band_end) * uFreq(band_end + 1)) ;

            inhBW_oct(dInt, cluster) = log2(inhBand(dInt, 2, cluster) / inhBand(dInt, 1, cluster));
            inhQ(dInt, cluster) = CF(cluster) / (inhBand(dInt, 2, cluster) - inhBand(dInt, 1, cluster));
        end
    end
    end
end

%% ----------------- Characterizing periodic effects ----------------- %%

%% OPTION 1: ISI characterization
% for vibrotactile and AM

% load relevant data
[cids, stimuli_parameters, aligned_spikes, ~, ~, ~, ~, StimResponseFiring] = loadData(OutPath, 3, BehaviorPath);

% inter spike interval analysis
[trials, clusters] = size(aligned_spikes);
aligned_spikes_ISI = cell(trials, clusters);
% set analysis window
if strcmp(stimuli_parameters.Par.Rec, "SxA") && strcmp(stimuli_parameters.Par.SomatosensoryWaveform, 'UniSine')
    % vibrotactile stimuli
    stimOn_ISIwindow = str2num(stimuli_parameters.Par.SomAudSOA)/1000;
    stimOff_ISIwindow = stimOn_ISIwindow + str2double(stimuli_parameters.Par.SomatosensoryStimTime)/1000;
elseif  strcmp(stimuli_parameters.Par.Rec, "AMn")
    % AM noise
    stimOn_ISIwindow = ((str2double(stimuli_parameters.Par.AMPostTime)) - (str2num(stimuli_parameters.Par.AMTransTime)))/1000;
    stimOff_ISIwindow = stimOn_ISIwindow + (str2num(stimuli_parameters.Par.AMTransTime))/1000;
end

for cluster = 1:clusters
    for trial = 1:trials
        % select stim on period only from aligned_spikes
        sel = aligned_spikes{trial, cluster} > stimOn_ISIwindow & aligned_spikes{trial, cluster} < stimOff_ISIwindow;
        taligned{trial, cluster} = aligned_spikes{trial,cluster}(sel);

        % log interspike interval
        aligned_spikes_ISI{trial, cluster} = log(diff(taligned{trial, cluster}));
    end
end

conditions = StimResponseFiring.conditions; % OO OA SO SA
uparamA = StimResponseFiring.amplitudes;
uparamB = StimResponseFiring.frequencies;
nparamA = length(uparamA);
nparamB = length(uparamB);
ISICVs = nan(nparamA, nparamB, length(conditions), length(cids)); % most freq ISI %amp,freq,con,cluster

spikeISI = cell(nparamA, nparamB, size(conditions, 2), length(cids));
for cluster = 1:length(StimResponseFiring.cids)
    for con = 1:length(conditions)
        for amp = 1:nparamA
            for freq = 1:nparamB

                % select unique stimulus combination
                if strcmp(stimuli_parameters.Par.Rec, "SxA")
                    index = stimuli_parameters.Stm.Amplitude == uparamA(amp) & stimuli_parameters.Stm.SomFreq == uparamB(freq) & strcmp(stimuli_parameters.Stm.MMType, conditions{con});
                    figTitle = [num2str(uparamB(freq)) 'Hz ' num2str(uparamA(amp)) 'V'];
                    figSGTitle = ['Unit: ' num2str(cids(cluster)) ', condition: ' conditions{con} '(' num2str(uparamA(amp)) 'V)'];
                elseif strcmp(stimuli_parameters.Par.Rec, "AMn")
                    index = stimuli_parameters.Stm.Intensity == uparamA(amp) & stimuli_parameters.Stm.Mf == uparamB(freq);
                    figTitle = [num2str(uparamB(freq)) 'Hz '];
                    figSGTitle = ['Unit: ' num2str(cids(cluster)) ', ' num2str(uparamA(amp)) 'dbSPL'];
                end

                if sum(index) == 0
                    continue
                end

                idx_rows = find(index);
                tspikeISI = [];

                % concatinate all ISIs
                for i = 1:sum(index)
                    ttISIspikes = aligned_spikes_ISI{idx_rows(i), cluster};
                    tspikeISI = [tspikeISI; ttISIspikes];
                end

                % to do: save SpikeISI in cell array to disentangle analysis from plotting
                spikeISI{amp, freq, con, cluster} = tspikeISI; %amp, freq, con, cluster

                % calculate coefficient of variation (CV = standard deviation / mean)
                ISICV = mean(exp(spikeISI{amp, freq, con, cluster}))/ std(exp(spikeISI{amp, freq, con, cluster}));
                ISICVs(amp,freq,con,cluster) = ISICV;

                % find most occuring ISI
                ISImax = max(N);
                ISImaxpos = find(N == ISImax);

                if length(ISImaxpos) == 1
                    ISIpeak = edges(ISImaxpos); % match with edges to get ISI
                else
                    %disp('two peaks found, check')
                    ISIpeak = edges(ISImaxpos(1)); % match with edges to get ISI
                end

                % save in matrix amp x freq x condition x cluster
                ISIpeaks(amp,freq,con,cluster) = ISIpeak;

            end
        end
    end
end

clear trial trials cluster clusters sel stimOn_ISIwindow stimOff_ISIwindow taligned tspikeISI ttspikeISI


%% OPTION 2: dataQuantification_TemporalCoding.m

session = 2;

[cids, stimuli_parameters, aligned_spikes, Srise, Sfall, sessions_TTLs, ~, ~, ~] = loadData(OutPath, session, BehaviorPath);


all_freqs = unique(stimuli_parameters.Stm.SomFreq);
amp = 0.3;
nClusters = length(cids);
cluster = find(cids == 441);

for freq = 1:length(all_freqs)
    
    index = strcmp(stimuli_parameters.Stm.MMType, "SO") & (stimuli_parameters.Stm.SomFreq == all_freqs(freq)) & (stimuli_parameters.Stm.Amplitude == amp);
    if sum(index) == 0
        continue
    end

    SOM_Hz = stimuli_parameters.Stm.SomFreq(index);

    Mf = unique(SOM_Hz);
    % ~~~ parameters for analysis ~~
    StimDur = 0.5; %s
    % StimDur = 1/Mf; % 1/Mf = 1 cycle
    SkipVal = 0;
    SkipMethod = '';
    nBins = 24;
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % calculation
    [CycT] = CycTimes(aligned_spikes(index, cluster),StimDur, Mf,SkipVal, SkipMethod);
    [Vs,Ph,Ray] = calcVS(CycT); % vector strength calculation

    % get Ph (phase) and plot

end
