%% Automated
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%Predefine

% Within folders, look for 'VR'
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        %Look for orientation tuning session
        sesopt = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'10*'));
        orisess=[];
        for sesidx=1:length(sesopt)
            filetmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},sesopt(sesidx).name,'*Parameters*.mat*'));
            if ~isempty(filetmp)
                tmpfile = load(fullfile(filetmp(1).folder,filetmp(1).name));
                if strcmp(tmpfile.parameters.Protocol.xfile,'stimGratingVariousParameters.x')
                    orisess=[orisess sesidx];
                end
            end
        end
                
       
        if isempty(orisess)
            disp(['No orientation tuning for ' fullfile(filetmp(1).folder)])
            continue
        end

        for sesidx=1:length(orisess)
            %% Loading parameters
            filetmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},sesopt(orisess(sesidx)).name,'*Parameters*.mat*'));
            tmpfile = load(fullfile(filetmp(1).folder,filetmp(1).name));
            pars = tmpfile.parameters.Protocol.pars;
            directionopt = pars(ismember(tmpfile.parameters.Protocol.parnames,'ori'),:);
            seq = tmpfile.parameters.Protocol.seqnums;
            seq = seq'; %1st all repetitions

            %Get direction per trial
            directiontrial = [];
            for rep = 1:tmpfile.parameters.Protocol.nrepeats
                directiontrial = [directiontrial directionopt(seq(rep,:)-(rep-1)*tmpfile.parameters.Protocol.npfilestimuli)];
            end
                        
            
            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(LocalDir,MiceOpt{midx},Dates4Mouse{didx});
            myLFDir = fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'ephys');
            sp = loadKSdir(myKsDir);
            
            % sp.st are spike times in seconds
            % sp.clu are cluster identities
            % spikes from clusters labeled "noise" have already been omitted
            
            %% basic quantification of spiking plot
            [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
            
            depthBins = 0:40:3840;
            ampBins = 0:30:min(max(spikeAmps),800);
            recordingDur = sp.st(end);
            
            [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
            
            
            %% Computing some useful details about spikes/neurons (like depths)
            
            [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
                templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
            
            %% synchronize data with timeline
            mySyncFile = dir(fullfile(myLFDir,'*','*nidq.bin'))
            expnr = strsplit(filetmp(1).folder,'\');
            expnr = expnr{end};
            [dataArray,samplerate,starttime,endtime] = ReadSGLXData(mySyncFile(1).name,mySyncFile(1).folder,expnr); %1 To sync with NP, 2 acquire live, 3 flipper
           
            
            % Sync Flippers
            tm = 1./nanmean(diff(AllTimeLines{midx,didx,orisess(sesidx)}(:,ismember(AllInputs,'timestamps')))); %Freq. Flipper
            Flippertimeline = AllTimeLines{midx,didx,orisess(sesidx)}(:,ismember(AllInputs,'flipper'));
            Flippertimeline(Flippertimeline<1.5)=0;
            Flippertimeline(Flippertimeline>1)=1;
            Actualtime = AllTimeLines{midx,didx,orisess(sesidx)}(:,ismember(AllInputs,'timestamps'));
                 
            FlipperGLX = dataArray(3,:);
            FlipperGLX = downsample(FlipperGLX,round(1./(tm./samplerate)));
            FlipperGLX(FlipperGLX<2)=0;
            FlipperGLX(FlipperGLX>2)=1;
            
            spikeTimestmp =  spikeTimes;
%             Clustr = sp.clu(spikeIndex);
%             spikeDepths = spikeDepths(spikeIndex);
            spikeTimesCorrected = nan(1,length(spikeTimestmp));
            figure;
            % align timelines every x seconds
            for i = 1:10:max(ceil(Actualtime))
                indx=find(Actualtime>=i & Actualtime<i+10);
                indx(indx>length(Actualtime))=[];
                indx(indx>length(FlipperGLX))=[];
                
                if isempty(indx)
                    break
                end
                
                if exist('h1')
                    delete(h1)
                    delete(h2)
                end
                % Check every x timepoints what the difference is between the
                % two flippers
                delayindx = finddelay(Flippertimeline(indx),FlipperGLX(indx));
                h1 = plot(Actualtime(indx),Flippertimeline(indx),'b-');
                hold on
                indx2 = indx+delayindx;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                h2 = plot(Actualtime(indx),FlipperGLX(indx2)+1,'r-');
                
                %delay in time
                delay = delayindx./tm;
                                
                drawnow
                indx = find(spikeTimestmp>=i+delay+starttime&spikeTimestmp<i+10+delay+starttime);
                if ~isempty(indx)
                    spikeTimesCorrected(indx)=spikeTimestmp(indx)-starttime-delay;
                end
                
            end
            
            %% Align to Trial Onset times
            tmp = AllTimeLines{midx,didx,orisess(sesidx)}(:,ismember(AllInputs,'photoDiode'));
            tmpmean = nanmean(tmp(:));
            
            tmp(tmp<tmpmean) = 0;
            tmp(tmp>=tmpmean)=1;
            tmpdiff = [0 abs(diff(tmp))'];
            tmpdiff = smooth(tmpdiff,round(tm));
            
            tmpdiff(tmpdiff>=quantile(tmpdiff(:),0.5))=1;
            tmpdiff(tmpdiff~=1)=0;
            
            
            starttrialidx = find(tmpdiff(1:end)~=0,1,'first'); %1st nonzero is when flipper starts
            endtrialidx = [];
            while starttrialidx(end) < length(tmpdiff)
                endtrialidx = [endtrialidx find(tmpdiff(starttrialidx(end):end)==0,1,'first')+starttrialidx(end)];
                ind = find(tmpdiff(endtrialidx(end):end)==1,1,'first');
                if isempty(ind)
                    break
                end
                starttrialidx = [starttrialidx ind+endtrialidx(end)];
                
            end
            
            ntrials = length(starttrialidx);
            SpikeTrialID = nan(1,length(spikeTimesCorrected));
            SpikeTrialTime = nan(1,length(spikeTimesCorrected));
            for trid = 1:ntrials
                tmp = find(spikeTimesCorrected>=Actualtime(starttrialidx(trid))&spikeTimesCorrected<=Actualtime(endtrialidx(trid)));
                SpikeTrialID(tmp)=trid;
                SpikeTrialTime(tmp) = spikeTimesCorrected(tmp)-Actualtime(starttrialidx(trid));
            end
            

            
            %% Looking at PSTHs aligned to some event
            
            % if you now have a vector of relevant event times, called eventTimes (but
            % not the cell array as above, just a vector):
            
            window = [-0.3 30]; % look at spike times from 0.3 sec before each event to 1 sec after
            
            % if your events come in different types, like different orientations of a
            % visual stimulus, then you can provide those values as "trial groups",
            % which will be used to construct a tuning curve. Here we just give a
            % vector of all ones.
            Events = [Actualtime(starttrialidx)];
            trialGroups = nan(1,length(Events));
            trialGroups(1:length(starttrialidx))=1;
            
            psthViewer(spikeTimesCorrected, Clustr, Events, window, trialGroups);
            
            % use left/right arrows to page through the clusters
            
            
            %% PSTHs across depth
            
            depthBinSize = 80; % in units of the channel coordinates, in this case µm
            timeBinSize = 0.01; % seconds
            bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
            psthType = 'norm'; % show the normalized version
            eventName = 'stimulus onset'; % for figure labeling
            
            [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimesCorrected, spikeDepths, ...
                depthBinSize, timeBinSize, Actualtime(starttrialidx), window, bslWin);
            
            figure;
            plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
            
            %% average to reward event
            % if you now have a vector of relevant event times, called eventTimes (but
            % not the cell array as above, just a vector):
            
            window = [-0.3 30]; % look at spike times from 0.3 sec before each event to 1 sec after
            
            % if your events come in different types, like different orientations of a
            % visual stimulus, then you can provide those values as "trial groups",
            % which will be used to construct a tuning curve. Here we just give a
            % vector of all ones.
            trialGroups = ones(size(Actualtime(starttrialidx)));
            
            psthViewer(spikeTimesCorrected, Clustr, Actualtime(starttrialidx), window, trialGroups);
            
            % use left/right arrows to page through the clusters
            
            
            %% Loading raw waveforms
            
            % To get the true waveforms of the spikes (not just kilosort's template
            % shapes), use the getWaveForms function:
            myKsRawDir = fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'ephys');
            apD = dir(fullfile(myKsRawDir,'*','*','*ap*.bin')); % AP band file from spikeGLX specifically
            
            gwfparams.dataDir = tmpfile(1).folder;    % KiloSort/Phy output folder
            gwfparams.fileName = apD(1).name;         % .dat file containing the raw
            gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
            gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
            gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = sp.clu(sp.clu==155);
            
            wf = getWaveForms(gwfparams);
            
            figure;
            imagesc(squeeze(wf.waveFormsMean))
            set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number');
            colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;

            
        end
        
    end
    
end