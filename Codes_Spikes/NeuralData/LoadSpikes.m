%% Automated
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
RedoAfterClustering=1;
Redo = 1; % Redo in general
%Predefine
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial
% Within folders, look for 'VR'
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    SessionOpttmp = SessionOpt(midx,:);
    SessionOpttmp(cellfun(@isempty,SessionOpttmp)) = [];
    for sesidx = 1:length(SessionOpttmp)
        thisdate = Dates4Mouse{find(~(cellfun(@isempty,cellfun(@(X) strfind(SessionOpttmp{sesidx},X),Dates4Mouse,'UniformOutput',0))))};
        thisses = strsplit(SessionOpttmp{sesidx},[thisdate '\']);
        thisses =thisses{end};
        if ~Redo && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'))
            if ~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'CuratedResults.mat'))
                continue
            elseif RedoAfterClustering
                myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
                myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                if isempty(myClusFile)
                    disp('This data is not yet curated with phy!!')
                    continue
                end
            end
        end
        
        if isempty(AllTimeLines{midx,sesidx})
            continue
        end
        
        %Saving directory
        if ~isdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses))
            mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses))
        end
        %% Loading data from kilosort/phy easily
        myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
        if ~isdir(myKsDir)
            continue
        end
        myLFDir = fullfile(DataDir,MiceOpt{midx},thisdate,'ephys');
        sp = loadKSdir(myKsDir);
        % sp.st are spike times in seconds
        % sp.clu are cluster identities
        % spikes from clusters labeled "noise" have already been omitted
        
        %% Plotting a driftmap
        [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
        figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'DriftMap.fig'))
        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'DriftMap.bmp'))
        %% basic quantification of spiking plot
        depthBins = 0:40:3840;
        ampBins = 0:30:min(max(spikeAmps),800);
        recordingDur = sp.st(end);
        
        [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
        plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);
        
        lfpFs = 2500;  % neuropixels phase3a - also 3b
        nChansInFile = 385;  % neuropixels phase3a, from spikeGLX
        %% Plotting some basics about LFPs
        if 0
            lfpD = dir(fullfile(myLFDir,'*','*','*.lf.bin')); % LFP file from spikeGLX specifically
            
            lfpFilename = fullfile(lfpD(1).folder, lfpD(1).name);
            
            
            [lfpByChannel, allPowerEst, F, allPowerVar] = ...
                lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);
            
            chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
            nC = length(chanMap);
            
            allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq
            
            % plot LFP power
            dispRange = [0 100]; % Hz
            marginalChans = [10:50:nC];
            freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};
            
            plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);
        end
        
        
        %% Computing some useful details about spikes/neurons (like depths)
        lfpD = dir(fullfile(myLFDir,'*','*','*.ap.bin')); % LFP file from spikeGLX specifically
        
        %             if length(lfpD)>1
        %                 lfpD = lfpD(find(~cellfun(@isempty,cellfun(@(X) strfind(X,'Concat'),{lfpD(:).folder},'UniformOutput',0))));
        %             end
        %             if length(lfpD)<1
        %                 keyboard
        %             end
        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
        spikeCluster = sp.clu;
        
        
           %% Get cluster information
        myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
        if isempty(myClusFile)
            disp('This data is not yet curated with phy!!')
            curratedflag=0;
            myClusFile = dir(fullfile(myKsDir,'cluster_group.tsv'));
            clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
            cluster_id = clusinfo.cluster_id;
            KSLabel = clusinfo.KSLabel;
            depth = nan(length(cluster_id),1);
            for clusid=1:length(depth)
                depth(clusid)=round(nanmean(spikeDepths(find(spikeCluster==clusid-1))));
            end
            myClusFile = dir(fullfile(myKsDir,'channel_map.npy'));
            channelmap = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
            myClusFile = dir(fullfile(myKsDir,'channel_positions.npy'));
            channelpos = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
            channelpos = channelpos(:,2);
            
            channel = nan(length(cluster_id),1);
            for clusid=1:length(channel)
                [minval,idx]=min(abs(depth(clusid)-channelpos));
                channel(clusid) = channelmap(idx);
            end
            clusinfo.ch = channel;
            clusinfo.depth = depth;
            clusinfo.id = cluster_id;
        else
            CurationDone = 1;
            save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'CuratedResults.mat'),'CurationDone')
            clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
            curratedflag=1;
            
            cluster_id = clusinfo.id;
            KSLabel = clusinfo.KSLabel;
            depth = clusinfo.depth;
            channel = clusinfo.ch;
        end
        Good_ID = find(sum(ismember(KSLabel,'good'),2)==4);
        
        %% Get Histology output
        histodone=0;
        if exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'),'dir')
            tmpfile = matfile(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'));
            try
                Depth2Area = tmpfile.Depth2Area;
                openfig(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'HistoEphysAlignment.fig'))
                histodone=1;
                disp('Data already aligned to histology')
            catch ME
            end
        end
        if ~histodone
            histofile = dir(fullfile(myKsDir,'channel_locations.json'));
            if isempty(histofile)
                histofile = dir(fullfile(myKsDir,'*.csv'));
                if length(histofile)>1
                    keyboard
                end
                if isempty(histofile)
                    
                    histofile = dir(fullfile(myKsDir,'probe_ccf.mat'));
                    
                    
                    if isempty(histofile)
                        disp('No channel locations known yet, do histology first')
                        histoflag=0;
                    else
                        
                        % alignment with petersen probe output
                        histoflag = 1;
                        histinfo =load(fullfile(histofile(1).folder,histofile(1).name));
                        
                          % Align ephys data with probe
                          [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,1,2);
                        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'HistoEphysAlignment.fig'))
                        
                    end
                else % Automatic alignment with brain globe output
                    histoflag = 1;
                    histinfo =readtable(fullfile(histofile(1).folder,histofile(1).name));
                    
                    % Align ephys data with probe
                    [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,0,1);

                    saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'HistoEphysAlignment.fig'))
                    
                end
            else
                histoflag=1;
                histinfo = fileread(fullfile(histofile(1).folder,histofile(1).name)); %Read json file
                histinfo = jsondecode(histinfo);% Decode json text
                
                
                % Align ephys data with probe
                [Depth2AreaPerChannel, Depth2AreaPerUnit]  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,0,1);
                saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'HistoEphysAlignment.fig'))
            end
        end
        
        
        
        %% load synchronization data
        % Timeline to go to:
        Actualtime = AllTimeLines{midx,sesidx}(:,ismember(AllInputs,'timestamps'));
        tmSR = 1./nanmean(diff(AllTimeLines{midx,sesidx}(:,ismember(AllInputs,'timestamps')))); %Freq. Flipper
        
        % Other events
        rewardevents = AllTimeLines{midx,sesidx}(:,ismember(AllInputs,'rewardEcho'));
        rewardevents(rewardevents<2)=0;
        rewardevents(rewardevents~=0)=1;
        RewardIndex = find(rewardevents);
        removevec = ones(1,length(rewardevents));
        removevec(RewardIndex) = [1; diff(RewardIndex)];
        RewardIndex = find(removevec>quantile(removevec(:),0.15));
        RewardTime = Actualtime(RewardIndex);
        
        % Sync Flippers
        Flippertimeline = AllTimeLines{midx,sesidx}(:,ismember(AllInputs,'flipper'));
        Flippertimeline(Flippertimeline<1.5)=0;
        Flippertimeline(Flippertimeline>1)=1;
        
        spikeTimestmp = spikeTimes;            % Spike times according to IMEC
        spikeTimesCorrected = nan(size(spikeTimestmp)); % corrected spike times (to timeline)
        warningflag = 0;
        % NIDQ FILE
        if nidq_sync_used(midx)
            mySyncFile = dir(fullfile(myLFDir,'*','*nidq.bin'))
            flag = 0;
            idx=1;
            while ~flag
                if idx>length(mySyncFile)
                    warning('No matching data, skip')
                    break
                end
                if ~isempty(strfind(mySyncFile(idx).name,'Concat'))
                    idx = idx+1;
                    continue
                end
                
                [dataArray,samplerate,starttime,endtime] = ReadSGLXData(mySyncFile(idx).name,mySyncFile(idx).folder,length(Actualtime)./tmSR); %1 To sync with NP, 2 acquire live, 3 flipper
                
                % Match this to the flipper signal
                FlipperGLX = dataArray(3,:);
                FlipperGLX = downsample(FlipperGLX,round(samplerate/tmSR)); %Make digital
                FlipperGLX(FlipperGLX<2)=0;
                FlipperGLX(FlipperGLX>2)=1;
                
                indx=1:length(Flippertimeline);
                % Initial estimate of delay
                delayestimate = finddelay(Flippertimeline(indx),FlipperGLX);
                indx2 = indx+delayestimate;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(FlipperGLX))=[];
                indx2(indx2>length(FlipperGLX))=[];
                if corr(Flippertimeline(indx),FlipperGLX(indx2)')<0.85 %Test if signals correlate enough
                    idx=idx+1;
                    continue
                end
                
                % NIDQ PULSE (=same clock as FlipperGLX)
                %                 mySyncFile = dir(fullfile(myLFDir,'*','*nidq.bin'));
                syncChanIndex = 1; %NPSync
                syncDatNIDQ = dataArray(syncChanIndex,:);%extractSyncChannel(mySyncFile(1).folder, 3, syncChanIndex); % Further digitize
                tmpmean = nanmean(syncDatNIDQ(:));
                syncDatNIDQ(syncDatNIDQ<tmpmean)=0;
                syncDatNIDQ(syncDatNIDQ>tmpmean)=1;
                %                 [nidqmeta] = (ReadMeta2(mySyncFile(1).folder));
                %Align to timeline
                syncDatNIDQ = downsample(syncDatNIDQ,round(samplerate/tmSR));
                %                                 syncDatNIDQTime = [1:length(syncDatNIDQ)]./str2double(nidqmeta.niSampRate);
                
                %IMEC PULSE
                syncChanIndex = 385;
                syncDatImec = extractSyncChannel(lfpD(idx).folder, nChansInFile, syncChanIndex);
                [Imecmeta] = (ReadMeta2(lfpD(idx).folder));
                
                % normalize between 0 and 1 (it's now in uint16 till 64)
                tmpmean = nanmean(syncDatImec(:));
                syncDatImec(syncDatImec<tmpmean)=0;
                syncDatImec(syncDatImec>tmpmean)=1;
                % Imec sync is not cut up in blocks yet; take correct
                %                 % sample
                %                 figure; plot(syncDatImec);
                %                 smoothsyncDatImec = smooth(single(syncDatImec));
                %                 hold on
                %                 startpoints = [];
                %                 endpoints = [];
                %                 tp = 1;
                %                 while tp<length(syncDatImec)
                %                     startpoints = [startpoints find(syncDatImec(tp:end)>0,1,'first')+tp-1];
                %                     line([startpoints(end) startpoints(end)],[0 1],'color',[0 1 0])
                %                     if isempty(find(syncDatImec(tp+1:end)>0,1,'first'))
                %                         break
                %                     end
                %                     tp = startpoints(end);
                %
                %                     endpoints = [endpoints find(syncDatImec(tp:end)<1,1,'first')+tp-1];
                %                     if isempty(endpoints)
                %                         endpoints = length(syncDataImec);
                %                     end
                %                     line([endpoints(end) endpoints(end)],[0 1],'color',[1 0 0])
                %                     if isempty(find(dataArray(2,tp+1:end)<1,1,'first'))
                %                         break
                %                     end
                %                     tp = endpoints(end);
                %
                %                 end
                %                 if length(endpoints)<length(startpoints)
                %                     endpoints = [endpoints size(dataArray,2)];
                %                     line([endpoints(end) endpoints(end)],[0 5],'color',[1 0 0])
                %
                %                 end
                %                 [~, SessionIndex] = min(abs(((endpoints-startpoints)./samplerate)-duration));
                %                 disp(['This is session ' num2str(SessionIndex)])
                %                 %Take two seconds extra if possible
                %                 startpoint = round(startpoints(SessionIndex)-samplerate* 10);
                %                 if startpoint<1
                %                     startpoint=1;
                %                 end
                %                 endpoint = round(endpoints(SessionIndex)+samplerate*10);
                %                 if endpoint>size(dataArray,2)
                %                     endpoint=size(dataArray,2);
                %                 end
                %                 dataArray = dataArray(:,startpoint:endpoint);
                %                 starttime = startpoint./ samplerate;
                %                 endtime = endpoint./ samplerate;
                
                
                
                %                 [dataArray,samplerate,starttime,endtime] = ReadSGLXData(mySyncFile(idx).name,mySyncFile(idx).folder,length(Actualtime)./tmSR); %1 To sync with NP, 2 acquire live, 3 flipper
                
                startidx = floor(starttime*str2double(Imecmeta.imSampRate));
                if startidx<1
                    startidx=1;
                end
                endidx = ceil(endtime*str2double(Imecmeta.imSampRate));
                if endidx>length(syncDatImec)
                    endidx=length(syncDatImec)
                end
                syncDatImec = syncDatImec(startidx:endidx);
                %Align to timeline
                syncDatImec = downsample(syncDatImec,round(str2double(Imecmeta.imSampRate)/tmSR));
                
                indx(indx2>length(syncDatImec))=[];
                indx2(indx2>length(syncDatImec))=[];
                % Initial estimate of delay
                delayestimate2 = finddelay(syncDatNIDQ(indx2),syncDatImec(indx2));
                indx3 = indx2+delayestimate2;
                indx(indx3<1) = [];
                indx2(indx3<1)=[];
                indx3(indx3<1)=[];
                indx(indx3>length(syncDatImec))=[];
                indx2(indx3>length(syncDatImec))=[];
                indx3(indx3>length(syncDatImec))=[];
                if corr(double(syncDatNIDQ(indx2)'),double(syncDatImec(indx3)'))<0.9
                    idx=idx+1;
                    continue
                end
                flag=1;
            end
            if ~flag
                continue
            end
            figure;
            per2check = 5;
            % align timelines every x seconds (align GLX to Timeline)
            for i = 1:per2check:max(ceil(Actualtime))
                indx=find(Actualtime>=i-1 & Actualtime<i+per2check-1);
                indx(indx>length(Actualtime))=[];
                indx(indx>length(FlipperGLX))=[];
                
                if isempty(indx)
                    break
                end
                
                if exist('h1')
                    delete(h1)
                    delete(h2)
                 
                end
                
                if exist('h3')
                       delete(h3)
                    delete(h4)
                end
                % Check every x timepoints what the difference is between the
                % two flippers
                indx2 = indx+delayestimate;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(FlipperGLX))=[];
                indx2(indx2>length(FlipperGLX))=[];
                if isempty(indx2)
                    continue
                end
                delayindx = finddelay(Flippertimeline(indx),FlipperGLX(indx2));
                
                %D = FINDDELAY(X,Y), where X and Y are row or column vectors of length
                %   LX and LY, respectively, returns an estimate of the delay D between X
                %   and Y, where X serves as the reference vector. If Y is delayed with
                %   respect to X, D is positive. If Y is advanced with respect to X, D is
                %   negative.
                
                indx2 = indx2+delayindx;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(FlipperGLX))=[];
                indx2(indx2>length(FlipperGLX))=[];
                if isempty(indx2)
                    continue
                end
                h1 = plot(Actualtime(indx),Flippertimeline(indx),'b-');
                hold on
                h2 = plot(Actualtime(indx),FlipperGLX(indx2)+1,'r-');
                
                
                indx3 = indx2+delayestimate2;
                indx(indx3<1)=[];
                indx2(indx3<1)=[];
                indx3(indx3<1)=[];
                indx(indx3>length(syncDatImec))=[];
                indx2(indx3>length(syncDatImec))=[];
                indx3(indx3>length(syncDatImec))=[];
                % sync channel data (ch. 385 NP) to Nidq channel 1
                delayindx2 =finddelay(syncDatNIDQ(indx2),syncDatImec(indx3));
                
                % +delayindx2
                indx3 = indx3+delayindx2;
                indx(indx3<1)=[];
                indx2(indx3<1)=[];
                indx3(indx3<1)=[];
                indx(indx3>length(syncDatImec))=[];
                indx2(indx3>length(syncDatImec))=[];
                indx3(indx3>length(syncDatImec))=[];
                
                h3 = plot(Actualtime(indx),syncDatNIDQ(indx2)+2,'r-');
                h4 = plot(Actualtime(indx),syncDatImec(indx3)+3,'k-');
                
                %sanity check
                %                     figure;
                %                     indx=abs((-delayindx2-delayestimate2-delayindx-delayestimate))+1:abs((-delayindx2-delayestimate2-delayindx-delayestimate))+5000;
                %                     plot(syncDatImec(indx),'k-')
                %                     hold on
                %                     plot(syncDatNIDQ(indx-delayindx2-delayestimate2)+1,'r-')
                %                     plot(FlipperGLX(indx-delayindx2-delayestimate2)+2,'r-')
                %                     plot(Flippertimeline(indx-delayindx2-delayestimate2-delayindx-delayestimate)+3,'b-')
                %
                %Convert back to time in imec data
                TL2ImecTime = (delayindx2+delayestimate2+delayindx+delayestimate)./tmSR;
                Imec2TLTime = (-delayindx2-delayestimate2-delayindx-delayestimate)./tmSR;
                
                %Find spikes in IMEC time, that fall in this window in
                %TIMELINE space
                spikeindx = find(spikeTimestmp>=TL2ImecTime+starttime+(i-1)&spikeTimestmp<=TL2ImecTime+starttime+(i+per2check-1));
                if ~isempty(spikeindx)
                    spikeTimesCorrected(spikeindx)=spikeTimestmp(spikeindx)-starttime+Imec2TLTime;
                end
                drawnow
            end
        else % Flipper directly in sync channel;
            
            indx=1:length(Flippertimeline);
            
            %IMEC PULSE
            syncChanIndex = 385;
            syncDatImec = extractSyncChannel(lfpD.folder, nChansInFile, syncChanIndex);
            [Imecmeta] = (ReadMeta2(lfpD.folder));
            
            % normalize between 0 and 1 (it's now in uint16 till 64)
            tmpmean = nanmean(syncDatImec(:));
            syncDatImec(syncDatImec<tmpmean)=0;
            syncDatImec(syncDatImec>tmpmean)=1;
            
            
            %Align to timeline
            syncDatImec = downsample(syncDatImec,round(str2double(Imecmeta.imSampRate)/tmSR));
            firstval = syncDatImec(1);
            minimdelay = find(syncDatImec~=firstval,1,'first');
            indx2 = indx+minimdelay;
            % Initial estimate of delay
            delayestimate = finddelay(Flippertimeline(indx),syncDatImec(indx2))+minimdelay;
            indx2 = indx+delayestimate;
            indx(indx2<1) = [];
            indx2(indx2<1)=[];
            indx(indx2>length(syncDatImec))=[];
            indx2(indx2>length(syncDatImec))=[];
            if corr(double(Flippertimeline(indx)),double(syncDatImec(indx2)'))<0.75
                warningflag = 1;
            end
            
            figure;
            per2check = 5;
            % align timelines every x seconds (align GLX to Timeline)
            for i = 1:per2check:max(ceil(Actualtime))
                indx=find(Actualtime>=i-1 & Actualtime<i+per2check-1);
                indx(indx>length(Actualtime))=[];
                indx(indx>length(syncDatImec))=[];
                
                if isempty(indx)
                    break
                end
                
                if exist('h1')
                    delete(h1)
                    delete(h2)
                end
                % Check every x timepoints what the difference is between the
                % two flippers
                indx2 = indx+delayestimate;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(syncDatImec))=[];
                indx2(indx2>length(syncDatImec))=[];
                if isempty(indx2)
                    continue
                end
                delayindx = finddelay(Flippertimeline(indx),syncDatImec(indx2));
                
                %D = FINDDELAY(X,Y), where X and Y are row or column vectors of length
                %   LX and LY, respectively, returns an estimate of the delay D between X
                %   and Y, where X serves as the reference vector. If Y is delayed with
                %   respect to X, D is positive. If Y is advanced with respect to X, D is
                %   negative.
                
                indx2 = indx2+delayindx;
                indx(indx2<1)=[];
                indx2(indx2<1)=[];
                indx(indx2>length(syncDatImec))=[];
                indx2(indx2>length(syncDatImec))=[];
                if isempty(indx2)
                    continue
                end
                h1 = plot(Actualtime(indx),Flippertimeline(indx),'b-');
                hold on
                h2 = plot(Actualtime(indx),syncDatImec(indx2)+1,'r-');
                
                
                %sanity check
                if warningflag
                    figure;
                    indx=abs((-delayindx-delayestimate))+1:abs((-delayindx-delayestimate))+5000;
                    plot(syncDatImec(indx),'k-')
                    hold on
                    plot(Flippertimeline(indx-delayindx-delayestimate)+1,'b-')
                    
                    AcceptableYN = input('Is this alignment acceptable? (y/n)','s');
                    while ~ismember(AcceptableYN,{'y','n'});
                        if strcmpi(AcceptableYN,'n')
                            disp('Not acceptable, skipping session')
                            continue
                        elseif strcmp(AcceptableYN,'y')
                            disp('Acceptable, continuing...')
                        end
                    end
                end
                if corr(double(Flippertimeline(indx)),double(syncDatImec(indx2)'))<0.9
                    keyboard
                end
                %Convert back to time in imec data
                TL2ImecTime = (delayindx+delayestimate)./tmSR;
                Imec2TLTime = (-delayindx-delayestimate)./tmSR;
                
                %Find spikes in IMEC time, that fall in this window in
                %TIMELINE space
                spikeindx = find(spikeTimestmp>=TL2ImecTime+(i-1)&spikeTimestmp<=TL2ImecTime+(i+per2check-1));
                if ~isempty(spikeindx)
                    spikeTimesCorrected(spikeindx)=spikeTimestmp(spikeindx)+Imec2TLTime;
                end
                drawnow
            end
        end
        
        if ~any(~isnan(spikeTimesCorrected))
            warning('No Spikes in this session... continue')
            continue
        end
        
        
        %% Align to Trial Onset times
        tmp = AllTimeLines{midx,sesidx}(:,ismember(AllInputs,'photoDiode'));
        tmp = (tmp-nanmin(tmp(:)))./(nanmax(tmp(:))-nanmin(tmp(:)));
        tresh=nanmedian(tmp);
        tmp(tmp>=tresh)=1;
        tmp(tmp<tresh) = 0;
        tmpdiff = [0 abs(diff(tmp))'];
        tresh = 0.05;
        while quantile(tmpdiff(:),tresh)==0
            tresh = tresh+0.01;
        end
        tmpdiff(tmpdiff>=quantile(tmpdiff(:),tresh))=1;
        tmpdiff(tmpdiff~=1)=0;
        
        
        starttrialidx = find(tmpdiff(1:end)~=0,1,'first'); %1st nonzero is when photodiode starts
        endtrialidx = [];
        trialid = 1;
        TrialDurations = nanmax(AllTime{midx,sesidx},[],2);
        while starttrialidx(end) < length(tmpdiff)
            %                 should be at least so many datapoints later:
            minendfr = find(Actualtime>Actualtime(starttrialidx(end)) + TrialDurations(trialid),1,'first');
            if isempty(minendfr)
                endtrialidx = [endtrialidx length(Actualtime)];
                break
            end
            endtrialidx = [endtrialidx find(tmpdiff(minendfr:end)==0,1,'first')+minendfr];
            ind = find(tmpdiff(endtrialidx(end):end)==1,1,'first');
            if isempty(ind) || ind+endtrialidx(end) == starttrialidx(end)
                break
            end
            if trialid+1>length(TrialDurations)
                warning('More onsets, but no more trials left?!')
                break
            end
            starttrialidx = [starttrialidx ind+endtrialidx(end)];
            trialid=trialid+1;
        end
        
        if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations>1/tmSR*10)
            warning('TrialDurations wrong, check this session')
            keyboard
        end
        
        ntrials = length(starttrialidx);
        
        
        SpikeTrialID = nan(1,length(spikeTimesCorrected));
        SpikeTrialTime = nan(1,length(spikeTimesCorrected));
        for trid = 1:ntrials
            if trid==1
                tmp = find(spikeTimesCorrected>=Actualtime(1)&spikeTimesCorrected<=Actualtime(starttrialidx(trid+1)));
            elseif trid==ntrials
                tmp = find(spikeTimesCorrected>=Actualtime(endtrialidx(trid-1))&spikeTimesCorrected<=Actualtime(end));
            else
                tmp = find(spikeTimesCorrected>=Actualtime(endtrialidx(trid-1))&spikeTimesCorrected<=Actualtime(starttrialidx(trid+1)));
            end
            SpikeTrialID(tmp)=trid;
            SpikeTrialTime(tmp) = spikeTimesCorrected(tmp)-Actualtime(starttrialidx(trid));
        end
        
        %remove NaNs
        spikeAmps(isnan(SpikeTrialID))=[];
        SpikeTrialTime(isnan(SpikeTrialID))=[];
        spikeCluster(isnan(SpikeTrialID)) = [];
        spikeDepths(isnan(SpikeTrialID))=[];
        spikeTimesCorrected(isnan(SpikeTrialID))=[];
        SpikeTrialID(isnan(SpikeTrialID))=[];
        
        %% Looking at PSTHs aligned to some event
        
        % if you now have a vector of relevant event times, called eventTimes (but
        % not the cell array as above, just a vector):
        
        window = [-2 30]; % look at spike times in this window
        
        if 0
            % if your events come in different types, like different orientations of a
            % visual stimulus, then you can provide those values as "trial groups",
            % which will be used to construct a tuning curve. Here we just give a
            % vector of all ones.
            Events = [Actualtime(starttrialidx);];
            trialGroups = nan(1,length(Events));
            trialGroups(1:length(starttrialidx))=1;
            
            psthViewer(spikeTimesCorrected, spikeCluster, Events, window, trialGroups);
            
            % use left/right arrows to page through the clusters
        end
        
        %% PSTHs across depth
        depthBinSize = 80; % in units of the channel coordinates, in this case µm
        timeBinSize = 0.01; % seconds
        bslWin = [-2 -0.05]; % window in which to compute "baseline" rates for normalization
        psthType = 'norm'; % show the normalized version
        eventName = 'stimulus onset'; % for figure labeling
        idx = ~isnan(spikeTimesCorrected);
        [timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimesCorrected(idx), spikeDepths(idx), ...
            depthBinSize, timeBinSize, Actualtime(starttrialidx), window, bslWin);
        
        figure;
        plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeRateDepth.fig'))
        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeRateDepth.bmp'))
        
        
     
        %% PSTH to position
        %             SpikeTrialTime(isnan(SpikeTrialID))=[];
        %             SpikeTrialID(isnan(SpikeTrialID))=[];
        disp('PSTH to Position')
        Pos = AllPosition{midx,sesidx};
        TrialTime = AllTime{midx,sesidx};
        timeBinSize = 0.05;
        bslWin = [-1,0];
        starttime = -2;
        endtime = 32;
        timevec = starttime+timeBinSize/2:timeBinSize:endtime-timeBinSize/2;
        %         SpikeRatePerPos = nan(ntrials,101,length(Good_ID));
        %         SpikeRatePerPosNorm = SpikeRatePerPos;
        SpikeRatePerTP = nan(ntrials,length(timevec),length(Good_ID));
        SpikeRatePerTPNorm = SpikeRatePerTP;
        for clusid = 1:length(Good_ID)
            clusid
            SpikesThisTrial = SpikeTrialTime(spikeCluster'== cluster_id(Good_ID(clusid)));
            if isempty(SpikesThisTrial)
                continue
            end
            
            for trid=1:ntrials
                SpikesThisTrial = SpikeTrialTime(SpikeTrialID == trid & spikeCluster'== cluster_id(Good_ID(clusid)));
                SpikesThisTrial=SpikesThisTrial(SpikesThisTrial>=starttime&SpikesThisTrial<=endtime);
                [N,edges]=histcounts(SpikesThisTrial,length(timevec));
                SpikeRatePerTP(trid,:,clusid) = N/timeBinSize; %spikes/second
                SpikeRatePerTPNorm(trid,:,clusid)=(SpikeRatePerTP(trid,:,clusid)-nanmean(SpikeRatePerTP(trid,timevec<0,clusid)))./nanmean(SpikeRatePerTP(trid,timevec<0,clusid));
                %                 parfor cmid=1:101
                %                     %find (first) time in trial where Pos was cmid
                %                     tmptime = TrialTime(trid,round(Pos(trid,:))>=-cmid+1 & round(Pos(trid,:))<-cmid+2);
                %                     if ~isempty(tmptime)
                %                         SpikeRatePerPos(trid,cmid,clusid)  = nanmean(SpikeRatePerTP(trid,timevec==min(tmptime),clusid));
                %                         SpikeRatePerPosNorm(trid,cmid,clusid)  = (nanmean(SpikeRatePerTP(trid,timevec==min(tmptime),clusid))-baseline)./baseline;
                %                     end
                %                 end
            end
        end
        
        depthtmp = depth(Good_ID);
        [sorteddepth,sortidx]=sort(depthtmp); % sort by depth
        
        %%
        % Show probe locations
        if 0
            figure;
            if histoflag
                subplot(2,5,[5,10])
                channeltmp = channel(Good_ID);
                channeltmp = channeltmp(sortidx)+1; %min channel = 0
                depthtmp = round(cell2mat(Depth2Area(1,channeltmp))); %Depth of units
                areatmp = Depth2Area(2,channeltmp); %Area of units
                coltmp = Depth2Area(3,channeltmp); % color of area for unit
                emptyidx = cellfun(@isempty,areatmp);
                areatmp(emptyidx)=[];
                depthtmp(emptyidx) = [];
                coltmp(emptyidx)=[];
                
                AreaOpt = unique(areatmp,'stable'); %identify areas
                for i=1:length(AreaOpt)
                    Idx = find(ismember(areatmp,AreaOpt{i}));
                    if strcmp(AreaOpt{i},'void')
                        h=patch([0 0 1 1],[Idx(1) Idx(end) Idx(end) Idx(1)],[1 1 1]);
                    else
                        h=patch([0 0 1 1],[Idx(1) Idx(end) Idx(end) Idx(1)],hex2rgb(coltmp{Idx(1)}{1}));
                    end
                    text(0.5,nanmean([Idx]),AreaOpt{i},'HorizontalAlignment','center')
                end
                set(gca,'ylim',[1,length(areatmp)],'YTickLabel',[])
            end
            %SpikeRate
            % plot average firing rate over position for different depths
            tmp = squeeze(nanmean(SpikeRatePerPosNorm(:,:,Good_ID),1));
            tmp = tmp(:,sortidx);            % in sorted depth
            subplot(2,5,[1:4,6:9])
            h=imagesc(tmp',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
            xlabel('Position (cm)')
            yticks = get(gca,'ytick');
            set(gca,'yTickLabel',sorteddepth(yticks))
            ylabel('Depth (mm)')
            colormap redblue
            set(h,'alphadata',~isnan(tmp)')
            box off
            colorbar
            title('SpikeRate (spikes/second)')
            saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeRatePositionAlongProbe.fig'))
            
            tmp = squeeze((SpikeRatePerPos(:,:,Good_ID)));
            tmp = tmp(:,:,sortidx);            % in sorted depth
            figure('name',[MiceOpt{midx} '\'  thisdate '\' thisses 'PSTH Per Area and Location'])
            for areaid=1:length(AreaOpt)
                subplot(ceil(sqrt(length(AreaOpt))),round(sqrt(length(AreaOpt))),areaid)
                idxopt = find(ismember(areatmp,AreaOpt{areaid}));
                hold on
                cols = distinguishable_colors(length(idxopt));
                for id = 1:length(idxopt)
                    shadedErrorBar(0:100,squeeze(nanmean(tmp(:,:,idxopt(id)),1)),...
                        squeeze(nanstd(tmp(:,:,idxopt(id)),[],1)./sqrt(size(tmp(:,:,idxopt(id)),1)-1)))
                end
                ylims = get(gca,'ylim');
                %             line([nanmean(AllRewardPos{midx,sesidx}) nanmean(AllRewardPos{midx,sesidx})],get(gca,'ylim'),'color',[0 1 0])
                patch([16 24 24 16],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
                patch([36 44 44 36],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
                patch([56 64 64 56],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
                patch([76 84 84 76],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
                
                title([AreaOpt{areaid}, ', ' num2str(sum(ismember(areatmp,AreaOpt{areaid}))) ' good units'])
            end
            suplabel('Position (cm)');
            suplabel('Average Firing Rate (spikes/sec)','y');
            
            saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'AvgSpikeRatePosition.fig'))
        end
        
        %% plot average firing rate over time
        figure;
        if histoflag
            channeltmp = channel(Good_ID);
            channeltmp = channeltmp(sortidx)+1; %min channel = 0
            depthtmp = round(cell2mat(Depth2Area(1,channeltmp))); %Depth of units
            areatmp = Depth2Area(2,channeltmp); %Area of units
            coltmp = Depth2Area(3,channeltmp); % color of area for unit
            emptyidx = cellfun(@isempty,areatmp);
            areatmp(emptyidx)=[];
            depthtmp(emptyidx) = [];
            coltmp(emptyidx)=[];
            
            %Order depth
            [depthtmp sortidx]=sort(depthtmp);
            areatmp = areatmp(sortidx);
            coltmp = coltmp(sortidx);
            
            [AreaOpt,IA,IC] = unique(areatmp,'stable'); %identify areas
            switchpoints = [1; find(diff(IC)~=0)+1; length(IC)]; %Find points where region changes
            AreaOpt = areatmp(switchpoints);
            subplot(2,5,[5,10])
            for i=1:length(AreaOpt)
                Idx = find(ismember(areatmp,AreaOpt{i}));
                if strcmp(AreaOpt{i},'void')
                    h=patch([0 0 1 1],[Idx(1) Idx(end) Idx(end) Idx(1)],[1 1 1]);
                else
                    h=patch([0 0 1 1],[Idx(1) Idx(end) Idx(end) Idx(1)],hex2rgb(coltmp{Idx(1)}{1}));
                end
                text(0.5,nanmean([Idx]),AreaOpt{i},'HorizontalAlignment','center')
            end
            set(gca,'ylim',[1,length(areatmp)],'YTickLabel',[])
        end
        %SpikeRate
        % plot average firing rate over position for different depths
        tmp = squeeze(nanmean(SpikeRatePerTP(:,:,sortidx),1));
        subplot(2,5,[1:4,6:9])
        h=imagesc(timevec,depthtmp,tmp',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
        xlabel('Time (seconds)')
        yticks = get(gca,'ytick');
        ylabel('Depth (mm)')
        colormap redblue
        set(h,'alphadata',~isnan(tmp)')
        box off
        colorbar
        title('Normalized SpikeRate (spikes/second)')
        set(gca,'ydir','normal')
        xlim([-1 31])
        saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeRateTimePOintAlongProbe.fig'))
        
%         
%         figure('name',[MiceOpt{midx} '\'  thisdate '\' thisses 'PSTH Per Area across Time'])
%         for areaid=1:length(AreaOpt)
%             subplot(ceil(sqrt(length(AreaOpt))),round(sqrt(length(AreaOpt))),areaid)
%             idxopt = find(ismember(areatmp,AreaOpt{areaid}));
%             hold on
%             cols = distinguishable_colors(length(idxopt));
%             for id = 1:length(idxopt)
%                 shadedErrorBar(timevec,squeeze(nanmean(SpikeRatePerTPNorm(:,:,idxopt(id)),1)),...
%                     squeeze(nanstd(SpikeRatePerTPNorm(:,:,idxopt(id)),[],1)./sqrt(size(SpikeRatePerTPNorm(:,:,idxopt(id)),1)-1)))
%             end
%             title([AreaOpt{areaid}, ', ' num2str(sum(ismember(areatmp,AreaOpt{areaid}))) ' good units'])
%             xlim([min(timevec) max(timevec)])
%         end
%         suplabel('time (sec)');
%         suplabel('Average Firing Rate (spikes/sec)','y');
%         
%         saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'AvgSpikeRateTime.fig'))
%         
        %% Save Preprocessed data
        save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'),'tmSR','SpikeTrialID','spikeCluster','SpikeTrialTime','SpikeRatePerTPNorm','SpikeRatePerTP','Depth2AreaPerChannel','Depth2AreaPerUnit','timevec','clusinfo','-v7.3')
        
        if 0
            %% Loading raw waveforms
            
            % To get the true waveforms of the spikes (not just kilosort's template
            % shapes), use the getWaveForms function:
            myKsRawDir = fullfile(DataDir,MiceOpt{midx},thisdate,'ephys');
            apD = dir(fullfile(myKsRawDir,'*','*','*ap*.bin')); % AP band file from spikeGLX specifically
            
            gwfparams.dataDir = apD(1).folder;    % KiloSort/Phy output folder
            gwfparams.dataDir2 = myKsDir;
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
        
        clear sp
        clear SpikeRatePerPos
        clear SpikeRatePerTP
        clear spikeTimesCorrected
        clear allP
        clear ampBins
        clear cdfs
        clear clusinfo
        clear dataArray
        clear FlipperGLX
        clear Flippertimeline
        clear removevec
        clear rewardevents
        clear spikeAmps
        clear spikeCluster
        clear spikeDepths
        clear spikeSites
        clear spikesThisTrial
        clear spikeTimes
        clear spikeTimestmp
        clear SpikeTrialID
        clear SpikeTrialTime
        clear tempsUnW
        clear Actualtime
        clear allPowerEst
        clear allPowerVar
        clear F
        close all
        clear starttime
        clear h1
        
    end
    
end