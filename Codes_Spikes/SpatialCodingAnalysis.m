%% Automated
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);


for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    SessionOpttmp = SessionOpt(midx,:);
    SessionOpttmp(cellfun(@isempty,SessionOpttmp)) = [];    
  
    %Initialize       
    AreaName = cell(0,0);
    uniquearean=0; %Number of unique areas in dataset    
    SpikeRatePerPos = cell(uniquearean,length(SessionOpttmp));
    SpikeRatePerTP = cell(uniquearean,length(SessionOpttmp));
    Timeline = cell(1,length(SessionOpttmp));
    
    for sesidx = 1:length(SessionOpttmp)
        thisdate = Dates4Mouse{find(~(cellfun(@isempty,cellfun(@(X) strfind(SessionOpttmp{sesidx},X),Dates4Mouse,'UniformOutput',0))))};
        thisses = strsplit(SessionOpttmp{sesidx},[thisdate '\']);
        thisses =thisses{end};
        if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'))
            disp(['Run LoadSpikes first for ' fullfile(SaveDir,MiceOpt{midx},thisdate,thisses)])
            continue
        else
            TMPDATA = load(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,'SpikeData.mat'));
        end        
        TrialTime = AllTime{midx,sesidx};

        %% Extract data -- good units only
        % Extract clusterinfo
        clusinfo = TMPDATA.clusinfo;
        if isfield(clusinfo,'id')
            cluster_id = clusinfo.id;            
        else
            cluster_id = clusinfo.cluster_id;          
        end
        try
        channel = clusinfo.ch;
        catch ME
            disp(ME)
            continue %for now
        end

        KSLabel = clusinfo.KSLabel;
        depth = clusinfo.depth;
        Good_ID = (sum(ismember(KSLabel,'good'),2)==4); %Identify good clusters
        
        % Extract depth to Area correspondance
        Depth2Area = TMPDATA.Depth2Area; % (1,:) = depth, (2,:) = area, (3,:) = colorcode for area
        
        channeltmp = channel(Good_ID)+1; % Find channel of good units
        depthtmp = round(cell2mat(Depth2Area(1,channeltmp))); %Depth of units
        areatmp = Depth2Area(2,channeltmp); %Area of units
        coltmp = Depth2Area(3,channeltmp); % color of area for unit
        emptyidx = cellfun(@isempty,areatmp);
        areatmp(emptyidx)={'Undefined'};
        depthtmp(emptyidx) =nan;
        coltmp(emptyidx)={{'#808080'}};
        
        Pos = AllPosition{midx,sesidx};
        TrialTime = AllTime{midx,sesidx};
        %% Plot Example cell/trial        
        SpikeTrialTime = TMPDATA.SpikeTrialTime;
        SpikeCluster = TMPDATA.spikeCluster;
        SpikeTrialID = TMPDATA.SpikeTrialID;
        timevec=TMPDATA.timevec;
        timesteps = min([nanmean((diff(timevec))) (diff(timevec))]);

        goodclusteridx = find(Good_ID);
        NrSpikesPerTP = nan(length(timevec),max(SpikeTrialID),length(goodclusteridx));
        NrSpikesPerPos = nan(101,max(SpikeTrialID),length(goodclusteridx));       
        
        %For frequency analysis
        Fs=100;
        L = length(linspace(-2,round(nanmax(abs(TrialTime(:)))),Fs));
        n=2^nextpow2(L);
        f = Fs*(0:(n/2))/n;
        PowerPerFrq= nan(length(f),length(goodclusteridx));
        
        cutofffreq = 8;
        for clusid = 1:length(goodclusteridx)
%             if exist('ThisClusterFig')
%                 delete(ThisClusterFig)
%             end
            
            ThisClusterFig = figure('name',['Cluster ' num2str(clusid)]);

            subid=1;
            subplot(3,4,subid)
            % get spikes & trial information
            thesespikes = SpikeTrialTime(SpikeCluster==goodclusteridx(clusid));
            thesetrials = SpikeTrialID(SpikeCluster==goodclusteridx(clusid));
            if isempty(thesespikes)
                continue
            end
            plot(thesespikes,thesetrials,'k.')
            xlim([-1 31])
            box off
            title(areatmp((clusid)))          
            ylabel('Trial')
            xlabel('Time (sec)')
            
            %Sort by trials
            spikespertrial = arrayfun(@(X) thesespikes(thesetrials==X),1:max(thesetrials),'UniformOutput',0);
            spikespertrial = cellfun(@(X) cell2mat(arrayfun(@(Y) sum(abs(X-Y)<timesteps),timevec,'UniformOutput',0)),spikespertrial,'UniformOutput',0);
            spikespertrial = cat(1,spikespertrial{:});
            spikespertrial(timevec>nanmax(TrialTime,[],2))=nan;            
            NrSpikesPerTP(:,:,clusid) = spikespertrial';
            
            if ~any(spikespertrial(:)~=0 & ~isnan(spikespertrial(:)))
                continue
            end
            % Average PSTHs
            subid = subid+1;
            subplot(3,4,subid)
            hold on     
            h(1) = plot(timevec,squeeze(nanmean(spikespertrial,1)),'k');
            xlim([-1 31])
            box off
            title([areatmp((clusid)) ' All'])
            xlabel('time (s)')
            ylabel('SpikeRate')           
            
            subid = subid+1;
            subplot(3,4,subid)
            hold on
             %Odd
            h(2) = plot(timevec,squeeze(nanmean(spikespertrial(1:2:end,:),1)),'r');
            xlim([-1 31])
            box off
            title([areatmp((clusid)) ' Odd'])
            xlabel('time (s)')
            ylabel('SpikeRate')
          
            subid = subid+1;
            subplot(3,4,subid)
            hold on
             %even
            h(3) = plot(timevec,squeeze(nanmean(spikespertrial(2:2:end,:),1)),'b');
            xlim([-1 31])
            box off
            title([areatmp((clusid)) ' even'])
            xlabel('time (s)')
            ylabel('SpikeRate')      
            
       
            %Position                    
            tmpidx = arrayfun(@(X) abs(Pos)>X-1&abs(Pos)<X,0:100,'UniformOutput',0); %find position spend there index
            tmpidx = cat(3,tmpidx{:}); %make matrix
            tmpPositionIndx = squeeze(num2cell(tmpidx,2));
            TimePerTrial = num2cell(TrialTime,2); %Make a cell for every trial out of Trialtime
            nrtrials = size(TimePerTrial,1);
            tmpspikesperposavg = nan(101,max(thesetrials));          
            distvec = 0:100;
            parfor trid = 1:nrtrials                
                timepointsperPos = cellfun(@(X) TimePerTrial{trid}(X),tmpPositionIndx(trid,:),'UniformOutput',0);
                idx = find(~cell2mat(cellfun(@isempty,timepointsperPos,'UniformOutput',0)));
                for cmid = 1:101 %non empty positions
                    %                     NrSpikesPerPos =
                    %                     nan(101,max(thesetrials),length(goodclusteridx));
                    %                     find spike index
                    if ismember(cmid,idx)
                        tmpspikesperposavg(cmid,trid) = nanmean(spikespertrial(trid,find(cell2mat(arrayfun(@(X) any(abs(timepointsperPos{cmid}-X)<timesteps),timevec,'UniformOutput',0)))));
                    end
                end     
            end          

            NrSpikesPerPos(:,:,clusid) = tmpspikesperposavg;
            
            subid = subid+1;
            subplot(3,4,subid)
            try
                h=imagesc(tmpspikesperposavg'.*-1,[quantile(tmpspikesperposavg(:),0.99).*-1 0]);
            catch
                h=imagesc(tmpspikesperposavg'.*-1);
            end
            set(h,'alphadata',~isnan(tmpspikesperposavg'))
            set(gca,'ydir','normal')
            colormap gray
            box off
            xlabel('Position (cm)')
            ylabel('trial')
            title(areatmp((clusid)))
                           
            subid = subid+1;
            subplot(3,4,subid)
            shadedErrorBar([0:100],squeeze(nanmean(tmpspikesperposavg,2)),nanstd(tmpspikesperposavg,[],2)./sqrt(sum(~isnan(tmpspikesperposavg),2)));
            box off
            title([areatmp((clusid)) ' all'])
            xlabel('Position (cm)')
            ylabel('SpikeRate')
            
            subid = subid+1;
            subplot(3,4,subid)
            %Odd
            shadedErrorBar([0:100],squeeze(nanmean(tmpspikesperposavg(:,1:2:end),2)),nanstd(tmpspikesperposavg(:,1:2:end),[],2)./sqrt(sum(~isnan(tmpspikesperposavg(:,1:2:end)),2)),'lineprops','r');
            box off
            title([areatmp((clusid)) ' Odd'])
            xlabel('Position (cm)')
            ylabel('SpikeRate')
            
            %even
            subid = subid+1;
            subplot(3,4,subid)
            shadedErrorBar([0:100],squeeze(nanmean(tmpspikesperposavg(:,2:2:end),2)),nanstd(tmpspikesperposavg(:,2:2:end),[],2)./sqrt(sum(~isnan(tmpspikesperposavg(:,2:2:end)),2)),'lineprops','b');
            box off
            title([areatmp((clusid)) ' Even'])
            xlabel('Position (cm)')
            ylabel('SpikeRate')
            
            % Frequency tagging?         
            %Convert data to nr spikes across time (every 1/100seconds), excluding spikes beyond a
            %trial
            tmptrials= thesetrials(thesespikes<=32&thesespikes>=-2);
            tmp = thesespikes(thesespikes<=32&thesespikes>-2); %tmpspikes
            psth = nan(length(linspace(-2,round(nanmax(abs(TrialTime(:)))),Fs)),max(tmptrials));
            for trid=1:max(tmptrials)
                psth(:,trid) = cell2mat(arrayfun(@(X) sum(abs(tmp(tmptrials==trid)-X)<1/Fs),linspace(-2,round(nanmax(abs(TrialTime(trid,:)))),Fs),'UniformOutput',0));                  
            end     
            avgpsth = nanmean(psth,2)./sqrt(sum(~isnan(psth),2));

            [wt,ftmp] = cwt(avgpsth,'Morse',Fs);            
            P = abs(wt).^2;
            relevantP = P(ftmp<cutofffreq,:);
            
            subid = subid+1;
            subplot(3,4,[subid subid+2])

            h=imagesc(linspace(-2,nanmax(abs(TrialTime(:))),Fs),ftmp(ftmp<cutofffreq),relevantP);
            xlabel('time (s)')
            ylabel('Frequency (Hz)')
            set(gca,'ydir','normal')
            box off
            title('TimeFrequency')
            colorbar
            
            subid = subid+3;
            Y=fft(avgpsth,n);
            P=abs(Y/n).^2;
            subplot(3,4,subid)
            h=plot(f,P(1:n/2+1));
            PowerPerFrq(:,clusid) = P(1:n/2+1);
            ylabel('|P(f)|^2')
            xlabel('Frequency (Hz)')
            xlim([2 cutofffreq])
            title('Power/Freq')

            drawnow
            
           if length(findobj('type','figure'))>100
             close all
           end
        end
      
        close all
        %% Save this per area
        % Add units to area specific cell
        areaopt = unique(areatmp);
        for areaid = 1:length(areaopt)
            if any(strcmp(AreaName,areaopt{areaid}))
                idx = find(strcmp(AreaName,areaopt{areaid}));
            else
                uniquearean = uniquearean+1; %another unique area
                idx = uniquearean;                
                AreaName{idx} = areaopt{areaid};
            end
            
            % Find all units in this recording in this area
            unitidx = find(ismember(areatmp,AreaName{idx}));
            SpikeRatePerTP{idx,sesidx} = NrSpikesPerTP(:,:,unitidx);
            
            Timeline{1,sesidx} = TMPDATA.timevec;
            timevec = TMPDATA.timevec;
            timesteps = [nanmean((diff(timevec))) (diff(timevec))];
          
            SpikeRatePerPos{idx,sesidx} =  NrSpikesPerPos(:,:,unitidx);
        end
    end
    save(fullfile(SaveDir,MiceOpt{midx},'SpatialCodingProcessed.mat'),'Timeline','AreaName','SpikeRatePerPos','SpikeRatePerTP')
    %% Z-score activity
    %     tmp = cellfun(@(X) X(:),SpikeRatePerPos,'UniformOutput',0);
    %     meanz = nanmean(cat(1,tmp{:}));
    %     sdz = nanstd(cat(1,tmp{:}));
    %     ZSRPerPos = cellfun(@(X) (X-meanz)./sdz,SpikeRatePerPos,'UniformOutput',0);
    %     ZSRPerPos = cellfun(@(X) (X-permute(repmat(nanmean(reshape(X,size(X,1)*size(X,2),[]),1)',[1,size(X,1),size(X,2)]),[2,3,1]))./...
    %         permute(repmat(nanstd(reshape(X,size(X,1)*size(X,2),[]),[],1)',[1,size(X,1),size(X,2)]),[2,3,1]),SpikeRatePerPosNorm,'UniformOutput',0);
    ZSRPerPos = cellfun(@(X) X,SpikeRatePerPos,'UniformOutput',0);
    
    ZSRPerPosOddTrials = cellfun(@(X) X(1:2:end,:,:),ZSRPerPos,'UniformOutput',0);
    ZSRPerPosEvenTrials = cellfun(@(X) X(2:2:end,:,:),ZSRPerPos,'UniformOutput',0);
    
    % Sort even trials on max firing rate (location)
    meanperunitEven = cellfun(@(X) squeeze(nanmean(X,1)),ZSRPerPosEvenTrials,'UniformOutput',0);
    meanperunitOdd = cellfun(@(X) squeeze(nanmean(X,1)),ZSRPerPosOddTrials,'UniformOutput',0);

    meanperunit = cellfun(@(X) squeeze(nanmean(X,1)),ZSRPerPos,'UniformOutput',0);
    semperunit = cellfun(@(X) squeeze(std(X,[],1)./sqrt(size(X,1)-1)),ZSRPerPos,'UniformOutput',0);

    %% For all unique areanames
    LineFig = figure('name','ZScorePerPositionPlots');
    ImgFigEven = figure('name','ZScorePerPositionEven');
    ImgFigOdd= figure('name','ZScorePerPositionOdd');

    for areaid = 1:uniquearean
        %Get average and sd
        try
            tmp = meanperunit(areaid,:);
            tmpmean = cat(2,tmp{~cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))});
            tmp = semperunit(areaid,:);
            tmpsem =  cat(2,tmp{~cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))});
        catch
            continue
        end
        if ~any(~isnan(tmpmean(:)))
            continue
        end
           %Plot fig
        figure(LineFig)
        subplot(ceil(sqrt(uniquearean)),round(sqrt(uniquearean)),areaid)
        hold on
        
        if size(tmpmean,1)==1
            tmpmean = tmpmean';
            tmpsem = tmpsem';
        end
%         for celid = 1:size(tmpmean,2)
%             shadedErrorBar(1:size(tmpmean,1),tmpmean(:,celid),tmpsem(:,celid),'lineProps',{'LineWidth',.5})
%         end
        shadedErrorBar(1:size(tmpmean,1),nanmean(tmpmean,2),nanstd(tmpmean,[],2)./sqrt(size(tmpmean,2)-1),'lineProps',{'color',[0 0 0],'LineWidth',1})
        title(AreaName{areaid})
        ylims = get(gca,'ylim');
%         line([nanmean(cat(2,AllRewardPos{midx,:})) nanmean(cat(2,AllRewardPos{midx,:}))],get(gca,'ylim'),'color',[0 1 0])
        
        patch([16 24 24 16],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([36 44 44 36],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([56 64 64 56],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([76 84 84 76],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        xlim([0 100])
        
        
        % EVEN FIGURE
        % Img Fig
        %Sort on maximum firing rate per neuron
        tmp = meanperunitEven(areaid,:);
        tmp(cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))) = [];
        tmp = cat(2,tmp{:});
        [~, maxid] = nanmax(tmp,[],1);
        [~,sortid] = sort(maxid); %Sort the neurons like this
        

        figure(ImgFigEven)
        subplot(ceil(sqrt(uniquearean)),round(sqrt(uniquearean)),areaid)
        h = imagesc(tmp(:,sortid)',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
        set(gca,'ydir','normal')
        colormap redblue
        set(h,'AlphaData',~isnan(tmpmean'))
        box off
        title(AreaName{areaid})
        
        % Odd Figure
        tmp = meanperunitOdd(areaid,:);
        tmp(cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))) = [];
        tmp = cat(2,tmp{:});
        
        figure(ImgFigOdd)
        subplot(ceil(sqrt(uniquearean)),round(sqrt(uniquearean)),areaid)
        h = imagesc(tmp(:,sortid)',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
        set(gca,'ydir','normal')
        colormap redblue
        set(h,'AlphaData',~isnan(tmpmean'))
        box off
        title(AreaName{areaid})
        
        
        
        
    end
        
    %% For predetermined areas
    LineFig = figure('name',[MiceOpt{midx} 'ZScorePerPositionPlots']);
    ImgFigEven = figure('name',[MiceOpt{midx} 'ZScorePerPositionEven']);
    ImgFigOdd= figure('name',[MiceOpt{midx} 'ZScorePerPositionOdd']);

    for areaid = 1:length(AREASOfInterest)
        %Get average and sd
        subareaidx = find(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,AREASOfInterest{areaid}),AreaName,'UniformOutput',0),'UniformOutput',0)));
        try
            tmp = meanperunit(subareaidx,:);
            tmp = tmp(:);            
            tmpmean = cat(2,tmp{~cell2mat(cellfun(@(X) size(X,1)<2,tmp,'UniformOutput',0))});
            tmp = semperunit(subareaidx,:);
            tmpsem =  cat(2,tmp{~cell2mat(cellfun(@(X) size(X,1)<2,tmp,'UniformOutput',0))});
        catch
            continue
        end
        if ~any(~isnan(tmpmean(:)))
            continue
        end
        %Plot fig
        figure(LineFig)
        subplot(ceil(sqrt(length(AREASOfInterest))),round(sqrt(length(AREASOfInterest))),areaid)
        hold on
        
        if size(tmpmean,1)==1
            tmpmean = tmpmean';
            tmpsem = tmpsem';
        end
        %         for celid = 1:size(tmpmean,2)
        %             shadedErrorBar(1:size(tmpmean,1),tmpmean(:,celid),tmpsem(:,celid),'lineProps',{'LineWidth',.5})
        %         end
        shadedErrorBar(1:size(tmpmean,1),nanmean(tmpmean,2),nanstd(tmpmean,[],2)./sqrt(size(tmpmean,2)-1),'lineProps',{'color',[0 0 0],'LineWidth',1})
        title([AREASOfInterest{areaid} ', n= ' num2str(size(tmpmean,2))])
        ylims = get(gca,'ylim');
%         line([nanmean(cat(2,AllRewardPos{midx,:})) nanmean(cat(2,AllRewardPos{midx,:}))],get(gca,'ylim'),'color',[0 1 0])
        
        patch([16 24 24 16],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([36 44 44 36],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([56 64 64 56],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        patch([76 84 84 76],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
        xlim([0 100])
        
        
        % EVEN FIGURE
        % Img Fig
        %Sort on maximum firing rate per neuron
        tmp = meanperunitEven(subareaidx,:);
        tmp = tmp(:);                    
        tmp(cell2mat(cellfun(@(X) size(X,1)<2,tmp,'UniformOutput',0))) = [];
        tmp = cat(2,tmp{:});
        [~, maxid] = nanmax(tmp,[],1);
        [~,sortid] = sort(maxid); %Sort the neurons like this
        
        
        figure(ImgFigEven)
        subplot(ceil(sqrt(length(AREASOfInterest))),round(sqrt(length(AREASOfInterest))),areaid)
        h = imagesc(tmp(:,sortid)',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
        set(gca,'ydir','normal')
        colormap redblue
        set(h,'AlphaData',~isnan(tmpmean'))
        box off
        title(AREASOfInterest{areaid})
        
        % Odd Figure
        tmp = meanperunitOdd(subareaidx,:);
        tmp = tmp(:);
        tmp(cell2mat(cellfun(@(X) size(X,1)<2,tmp,'UniformOutput',0))) = [];
        tmp = cat(2,tmp{:});
        
        figure(ImgFigOdd)
        subplot(ceil(sqrt(length(AREASOfInterest))),round(sqrt(length(AREASOfInterest))),areaid)
        h = imagesc(tmp(:,sortid)',[-quantile(abs(tmp(:)),0.95) quantile(abs(tmp(:)),0.95)]);
        set(gca,'ydir','normal')
        colormap redblue
        set(h,'AlphaData',~isnan(tmpmean'))
        box off
        title(AREASOfInterest{areaid})
  
    end
    
    %% Frequency tagging?
    try
    Fs = TMPDATA.tmSR; %Sampling Freq
    catch
        Fs=1000;
    end
    T=1/Fs; % Sampling period    
    
    
    L = max(TMPDATA.SpikeTrialTime)*Fs;
    t = (0:L-1)*T; %time vector
    n=2^nextpow2(L);
    
    tic
    Data = arrayfun(@(X) sum((TMPDATA.SpikeTrialTime-X)<T),t,'UniformOutput',0);
    toc
    Data = cell2mat(Data);
    Y = fft(Data,n);
    P = abs(Y/n).^2;
    
    figure;
    plot(Fs*(0:(n/2))/n,log(P(1:n/2+1)));
    xlabel('Frequency (f)')
    ylabel('|P(f)|^2')
    xlim([1 11])
    
    
    % Z-score activity
    %     tmp = cellfun(@(X) X(:),SpikeRatePerPos,'UniformOutput',0);
    %     meanz = nanmean(cat(1,tmp{:}));
    %     sdz = nanstd(cat(1,tmp{:}));
    %     ZSRPerPos = cellfun(@(X) (X-meanz)./sdz,SpikeRatePerPos,'UniformOutput',0);
%     SpikeRatePerTPNorm =  cell(uniquearean,length(SessionOpttmp));
    timevec = unique([Timeline{:}]);%unique timeline (should be same timeline for all sessions!)
    T = unique(round(diff(timevec)*100)/100); %sampling frequency
    Fs = 1./T; %Sampling period
    f = Fs.*(0:(length(timevec)/2))./length(timevec);

    fftY = cellfun(@(X) fft(reshape(permute(X,[2,1,3]),size(X,1)*size(X,2),[]),[],1),SpikeRatePerTPNorm,'UniformOutput',0); %fast fourier transform
    Power = cellfun(@(X) abs(X/length(timevec)).^2,fftY,'UniformOutput',0); %Power

    
    %Plot(fftY,P1)
    %title('Single-Sided Amplitude Spectrum of S(t)')
    %xlabel('f (Hz)')
    %ylabel('|P1(f)|')
    
    LineFig = figure('name',[MiceOpt{midx} 'FrequencyTagging']);   
    for areaid = 1:length(AREASOfInterest)
        %Get average and sd
        subareaidx = find(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,AREASOfInterest{areaid}),AreaName,'UniformOutput',0),'UniformOutput',0)));
        try
            tmp = Power(subareaidx,:);
            tmp = tmp(:);            
            tmp(cell2mat(cellfun(@(X) size(X,1)<2,tmp,'UniformOutput',0)))=[];
            tmp=cellfun(@(X) X(1:length(timevec)/2+1,:),tmp,'UniformOutput',0);
            tmpmean = cat(2,tmp{:});
   
        catch
            continue
        end
        if ~any(~isnan(tmpmean(:)))
            continue
        end
        %Plot fig
        figure(LineFig)
        subplot(ceil(sqrt(length(AREASOfInterest))),round(sqrt(length(AREASOfInterest))),areaid)
        hold on
        
        if size(tmpmean,1)==1
            tmpmean = tmpmean';
        end
        %         for celid = 1:size(tmpmean,2)
        %             shadedErrorBar(1:size(tmpmean,1),tmpmean(:,celid),tmpsem(:,celid),'lineProps',{'LineWidth',.5})
        %         end
        shadedErrorBar(f,nanmean(tmpmean,2),nanstd(tmpmean,[],2)./sqrt(size(tmpmean,2)-1),'lineProps',{'color',[0 0 0],'LineWidth',1})
        title([AREASOfInterest{areaid} ', n= ' num2str(size(tmpmean,2))])
        xlim([1 10]);
        %         line([nanmean(cat(2,AllRewardPos{midx,:})) nanmean(cat(2,AllRewardPos{midx,:}))],get(gca,'ylim'),'color',[0 1 0])
        xlabel('Frequency (Hz)')
        ylabel('|P(f)|^2')
        
        
    end
    
end

%% 