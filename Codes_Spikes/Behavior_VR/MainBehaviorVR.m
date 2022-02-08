function [AllContr,AllSound,AllPosition,AllTime,AllTraj,AllRewarded,AllRewardPos,SessionOpt] = MainBehaviorVR(DataDir,MiceOpt,maxsessnr,MinDist2Include)


%% Behavior of VR

%% Automated
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%Predefine
SessionOpt = cell(length(MiceOpt),max(cell2mat(cellfun(@length,DateOpt,'UniformOutput',0)))+maxsessnr);
AllContr = cell(length(MiceOpt),max(cell2mat(cellfun(@length,DateOpt,'UniformOutput',0))));
AllSound = AllContr;
AllPosition = AllContr;
AllTime = AllContr;
AllTraj = AllContr;
AllRewarded = AllContr;
AllRewardPos = AllContr;
nrsess=nan(1,length(MiceOpt));
% Within folders, look for 'VR'
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    sescount = 1;
    for didx = 1:length(Dates4Mouse)
        SessionOpttmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'*10*'));%All my sessions start with 10(1/2/3/etc.)
        for sesidx=1:length(SessionOpttmp)
            %Predefine
            Contr = [];
            Sound = [];
            Position = [];
            Time = [];
            Traj = [];
            Rewarded = [];
            RewardPos = [];
            flag = 0;
            VRDataFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*VRBehavior*'));
            
            if isempty(VRDataFiles)
                continue %Not VR, continue
            else
                flag = 1;
            end
            
            %Store
            SessionOpt(midx,sescount) = {fullfile(SessionOpttmp(sesidx).folder,SessionOpttmp(sesidx).name)};            
            % Latest trial has all information you need
            for tridx = length(VRDataFiles)
                tmp = load(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,VRDataFiles(tridx).name));
                
                %Save variables
                nrtrials = tmp.TRIAL.info.no;
                Contr = [Contr tmp.TRIAL.trialContr(1:nrtrials)];
                Sound = [Sound tmp.TRIAL.Sound(1:nrtrials)];
                Positiontmp = tmp.TRIAL.posdata(1:nrtrials,:,3);
                Positiontmp(Positiontmp==0)=nan;
                Position = cat(1,Position, Positiontmp);
                Timetmp = tmp.TRIAL.time(1:nrtrials,:);
                Timetmp(Timetmp==0)=nan;
                Timetmp = Timetmp - Timetmp(:,2); %Relative to trial onset
                Time = cat(1,Time, Timetmp);
                Trajtmp = tmp.TRIAL.traj(1:nrtrials,:);
                Trajtmp(Trajtmp==0)=nan;
                Traj = cat(1,Traj, Trajtmp);
                Rewardedtmp = zeros(1,nrtrials);
                Rewardedtmp(tmp.REWARD.TRIAL)=1;
                Rewarded = [Rewarded Rewardedtmp];
                RewardPos = [RewardPos tmp.TRIAL.trialRewPos];
            end
            if flag
                AllContr{midx,sescount} = Contr;
                AllSound{midx,sescount} = Sound;
                AllPosition{midx,sescount} = Position;
                AllTime{midx,sescount} = Time;
                AllTraj{midx,sescount} = Traj;
                AllRewarded{midx,sescount}= Rewarded;
                AllRewardPos{midx,sescount}=RewardPos;
                sescount = sescount+1;
            end
        end
        
    end
    nrsess(midx) = sescount;
    
end


%% Trial vs. Distance
figure('name','Trial vs Distance');
for midx=1:length(MiceOpt)
      dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    
    tmpdist = cat(1,AllTraj{midx,dayidx});
    tmptime = cat(1,AllTime{midx,dayidx});
    tmprewpos = cat(2,AllRewardPos{midx,dayidx});
    tmprew = cat(2,AllRewarded{midx,dayidx});
    trid = repmat(1:size(tmpdist,1),size(tmpdist,2),1)';
    for tidx=1:size(trid,1)
        scatter(tmpdist(tidx,:),trid(tidx,:),14,tmptime(tidx,:))
        hold on
        if tmprew(tidx) %Draw in reward
            line([tmprewpos(tidx) tmprewpos(tidx)],[tidx-0.5 tidx+0.5],'color',[0 1 0],'LineWidth',3)
        end
    end
    xlabel('Distance (cm)')
    ylabel('Trial')
    hold on
    
    title(MiceOpt{midx})
    h = colorbar;
    set(get(h,'label'),'string','time (s)');
    colormap hot
    xlim([0 100])
    ylim([1 max(trid(:))])
    
end

%% Look at reward time / problem?
figure('name','Reward Time refresh differences');
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    
    subplot(length(MiceOpt),1,midx)
    tmpdist = cat(1,AllTraj{midx,dayidx});
    tmptime = cat(1,AllTime{midx,dayidx});
    tmprewardpos = cat(2,AllRewardPos{midx,dayidx});
    tmpreward = cat(2,AllRewarded{midx,dayidx});
    trid = repmat(1:size(tmpdist,1),size(tmpdist,2),1)';
    
    diffreward = [];
    diffnormal = [];
    for tidx=1:size(trid,1)
        %find last timepoint point before 70
        if tmpreward(tidx)
            afterrewardidx = find(tmpdist(tidx,:)>tmprewardpos(tidx),1,'first');
        else
            continue
        end
        if isempty(afterrewardidx)
            continue
        end
        beforerewardidx=find(tmpdist(tidx,:)<tmprewardpos(tidx),1,'last');
        if (tmptime(tidx,afterrewardidx)-tmptime(tidx,beforerewardidx))<0
            continue
        end
        diffreward =[diffreward (tmptime(tidx,afterrewardidx)-tmptime(tidx,beforerewardidx))];
        diffnormal = [diffnormal [diff(tmptime(tidx,:))]];
    end
    hist(diffnormal(~isnan(diffnormal)),1000)
    hold on
    line([quantile(diffnormal(~isnan(diffnormal)),0.95) quantile(diffnormal(~isnan(diffnormal)),0.95)],get(gca,'ylim'))
    line([nanmean(diffreward) nanmean(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','-')
    line([nanmean(diffreward)-nanstd(diffreward) nanmean(diffreward)-nanstd(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    line([nanmean(diffreward)+nanstd(diffreward) nanmean(diffreward)+nanstd(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    
    xlim([quantile(diffnormal(~isnan(diffnormal)),0.005) quantile(diffnormal(~isnan(diffnormal)),0.995)])
    
    xlabel('time difference')
    ylabel('nr. datapoints')
    hold on
    
    title(MiceOpt{midx})
end

figure('name','refresh difference larger than 0.032');
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    
    subplot(length(MiceOpt),1,midx)
    tmpdist = cat(1,AllTraj{midx,dayidx});
    tmptime = cat(1,AllTime{midx,dayidx});
    trid = repmat(1:size(tmpdist,1),size(tmpdist,2),1)';
    difftime = [];
    for tidx=1:size(trid,1)
        %find timepoint where timediff is larger than 0.032
        
        difftime = [difftime tmptime(tidx,find([nan [diff(tmpdist(tidx,:))]]>0.032))];
        
    end
    hist(difftime(~isnan(difftime)),1000)
    hold on
    line([quantile(difftime(~isnan(difftime)),0.95) quantile(difftime(~isnan(difftime)),0.95)],get(gca,'ylim'))
    line([nanmean(diffreward) nanmean(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','-')
    line([nanmean(diffreward)-nanstd(diffreward) nanmean(diffreward)-nanstd(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    line([nanmean(diffreward)+nanstd(diffreward) nanmean(diffreward)+nanstd(diffreward)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    
    xlim([quantile(diffnormal(~isnan(diffnormal)),0.01) quantile(diffnormal(~isnan(diffnormal)),0.99)])
    
    xlabel('time difference')
    ylabel('nr. datapoints')
    hold on
    
    title(MiceOpt{midx})
end

%% Distance Travelled across contrasts
ContrOpt = unique([AllContr{:}]);
AudiOpt = unique([AllSound{:}]);
meancontr = nan(length(MiceOpt),length(ContrOpt),length(AudiOpt));
semcontr = meancontr;
legendtext = {};
for aidx=1:length(AudiOpt)
    legendtext = {legendtext{:} ['Auditory = ' num2str(AudiOpt(aidx))]};
end

figure('name','DistanceWithContrast')
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    subplot(length(MiceOpt),1,midx)
    
    DistTot = max(cat(1,AllTraj{midx,dayidx}),[],2);
    Contr = cat(2,AllContr{midx,dayidx})';
    Audi = cat(2,AllSound{midx,dayidx})';
    for cidx=1:length(ContrOpt)
        for aidx=1:length(AudiOpt)
            meancontr(midx,cidx,aidx) = nanmean(DistTot(Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx)));
            semcontr(midx,cidx,aidx)= nanstd(DistTot(Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx)))./sqrt(sum(Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx))-1);
        end
    end
    barwitherr(squeeze(semcontr(midx,:,:)),squeeze(meancontr(midx,:,:)))
    
    box off
    set(gca,'XTickLabel',ContrOpt)
    title(MiceOpt{midx})
    ylabel('Distance')
    xlabel('Contrast')
    legend(legendtext)
    
end

%% Same across mice
figure('name','DistanceWithContrastAvg')
h = barwitherr(squeeze(nanstd(meancontr,[],1)./sqrt(length(MiceOpt)-1)),squeeze(nanmean(meancontr,1)));
h(1).FaceColor = [0 0 0];
h(2).FaceColor=[0.5 0.5 0.5];
legend(legendtext)

box off
set(gca,'XTickLabel',ContrOpt)
title('Across mice')
ylabel('Distance (cm)')
xlabel('Contrast')

makepretty
g_contr = repmat(ContrOpt',[1,length(AudiOpt),length(MiceOpt)]);
g_contr = permute(g_contr,[3,1,2]);
g_audi = repmat(AudiOpt',[1,length(ContrOpt),length(MiceOpt)]);
g_audi = permute(g_audi,[3,2,1]);
p= anovan(meancontr(:),{g_contr(:),g_audi(:)},'model','interaction','varnames',{'Visual Contrast','Auditory'});



%% Trial Selection
TrialSelection = cellfun(@(X) nanmax(X,[],2)>MinDist2Include, AllTraj,'UniformOutput',0);

%% Speed acrosss contrasts
TimeDiff = cellfun(@(X) cat(2,nan(size(X,1),1),diff(X,[],2)),AllTime,'UniformOutput',0);
PosDiff = cellfun(@(X) cat(2,nan(size(X,1),1),diff(X,[],2)),AllTraj,'UniformOutput',0);

Speed = cellfun(@(X,Y) Y./X,TimeDiff,PosDiff,'UniformOutput',0);
meancontr = nan(length(MiceOpt),length(ContrOpt),length(AudiOpt));
semcontr = meancontr;
legendtext = {};
for aidx=1:length(AudiOpt)
    legendtext = {legendtext{:} ['Auditory = ' num2str(AudiOpt(aidx))]};
end

figure('name','SpeedWithContrast')
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    
    DistTot = max(cat(1,Speed{midx,dayidx}),[],2);
    Contr = cat(2,AllContr{midx,dayidx})';
    Audi = cat(2,AllSound{midx,dayidx})';
    Trials2Include = cat(1,TrialSelection{midx,dayidx});
    for cidx=1:length(ContrOpt)
        for aidx=1:length(AudiOpt)
            meancontr(midx,cidx,aidx) = nanmean(DistTot(Trials2Include&Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx)));
            semcontr(midx,cidx,aidx)= nanstd(DistTot(Trials2Include&Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx)))./sqrt(sum(Trials2Include&Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx))-1);
        end
    end
    barwitherr(squeeze(semcontr(midx,:,:)),squeeze(meancontr(midx,:,:)))
    
    box off
    set(gca,'XTickLabel',ContrOpt)
    title(MiceOpt{midx})
    ylabel('Speed (cm/sec)')
    xlabel('Contrast')
    makepretty
end
legend(legendtext)
makepretty

%% Same across mice
figure('name','SpeedWithContrastAvg')
h = barwitherr(squeeze(nanstd(meancontr,[],1)./sqrt(length(MiceOpt)-1)),squeeze(nanmean(meancontr,1)));
h(1).FaceColor = [0 0 0];
h(2).FaceColor=[0.5 0.5 0.5];
legend(legendtext)

box off
set(gca,'XTickLabel',ContrOpt)
title('Across mice')
ylabel('Speed (cm/sec)')
xlabel('Contrast')

makepretty
g_contr = repmat(ContrOpt',[1,length(AudiOpt),length(MiceOpt)]);
g_contr = permute(g_contr,[3,1,2]);
g_audi = repmat(AudiOpt',[1,length(ContrOpt),length(MiceOpt)]);
g_audi = permute(g_audi,[3,2,1]);
p= anovan(meancontr(:),{g_contr(:),g_audi(:)},'model','interaction','varnames',{'Visual Contrast','Auditory'});

%% Speed vs location
cols = lines(length(MiceOpt));
figure('name','SpeedWithLocation')
clear h
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    
    Position = cat(1,AllTraj{midx,dayidx});
    Speedtmp = cat(1,Speed{midx,dayidx});
    Trials2Include = cat(1,TrialSelection{midx,dayidx});
    Position = Position(Trials2Include,:);
    Speedtmp=Speedtmp(Trials2Include,:);
    
    PositionVec = 0:1:100;
    SpeedsPerLoc = nan(length(PositionVec),2); %Mean, Sem
    %Find speeds at these position
    for posid=1:length(PositionVec)
        tmpspeed = (Speedtmp(round(Position)==PositionVec(posid)));
        if isempty(tmpspeed)
            continue
        end
        SpeedsPerLoc(posid,1)=nanmean(tmpspeed);
        SpeedsPerLoc(posid,2)=nanstd(tmpspeed)./sqrt(length(tmpspeed)-1);
    end
    
    h(midx) = shadedErrorBar(PositionVec,SpeedsPerLoc(:,1),SpeedsPerLoc(:,2),'lineprops',{'color',cols(midx,:)});
    hold on
    
end
ylims = get(gca,'ylim');
line([70 70],get(gca,'ylim'),'color',[0 1 0])
patch([16 24 24 16],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
patch([36 44 44 36],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
patch([56 64 64 56],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
patch([76 84 84 76],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
ylim(ylims)

box off
ylabel('Speed (cm/sec)')
xlabel('Position (cm)')
legend([h(:).mainLine],MiceOpt)
makepretty

%% Same, split up by condition
cols = lines(length(ContrOpt));
Linestls = {'-','--'};
figure('name','SpeedWithLocationAcross Conditions')
for midx=1:length(MiceOpt)
    dayidx=find(~cell2mat(cellfun(@isempty,AllContr(midx,:),'UniformOutput',0)));
    
    
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    
    DistTot = max(cat(1,AllTraj{midx,dayidx}),[],2);
    Contr = cat(2,AllContr{midx,dayidx})';
    Audi = cat(2,AllSound{midx,dayidx})';
    Position = cat(1,AllTraj{midx,dayidx});
    Speedtmp = cat(1,Speed{midx,dayidx});
    Trials2Include = cat(1,TrialSelection{midx,dayidx});
    
    DistTot = DistTot(Trials2Include);
    Contr = Contr(Trials2Include);
    Audi=Audi(Trials2Include);
    
    Position = Position(Trials2Include,:);
    Speedtmp=Speedtmp(Trials2Include,:);
    
    PositionVec = 0:1:100;
    SpeedsPerLoc = nan(length(PositionVec),2,length(ContrOpt),length(AudiOpt)); %Mean, Sem
    %Find speeds at these position
    for cidx=1:length(ContrOpt)
        for aidx=1:length(AudiOpt)
            for posid=1:length(PositionVec)
                tmpspeed = Speedtmp(round(Position)==PositionVec(posid)&Contr==ContrOpt(cidx)&Audi==AudiOpt(aidx));
                if isempty(tmpspeed)
                    continue
                end
                SpeedsPerLoc(posid,1,cidx,aidx)=nanmean(tmpspeed);
                SpeedsPerLoc(posid,2,cidx,aidx)=nanstd(tmpspeed)./sqrt(length(tmpspeed)-1);
            end
            h(cidx) = shadedErrorBar(PositionVec,SpeedsPerLoc(:,1,cidx,aidx),SpeedsPerLoc(:,2,cidx,aidx),'lineprops',{'color',cols(cidx,:),'LineStyle',Linestls{aidx}});
            hold on
        end
    end
    
    
    
    ylims = get(gca,'ylim');
    %     line([70 70],get(gca,'ylim'),'color',[0 1 0])
    patch([16 24 24 16],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
    patch([36 44 44 36],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
    patch([56 64 64 56],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
    patch([76 84 84 76],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
    title(MiceOpt{midx})
    box off
    ylabel('Speed (cm/sec)')
    xlabel('Position (cm)')
    Legendname = {};
    for cidx=1:length(ContrOpt)
        Legendname = {Legendname{:} ['c=' num2str(ContrOpt(cidx))]};
    end
    makepretty

end
legend([h(:).mainLine],Legendname)
makepretty


