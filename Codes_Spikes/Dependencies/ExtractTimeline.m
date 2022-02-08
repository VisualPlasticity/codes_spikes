function [AllTimeLines,AllInputs] = ExtractTimeline(DataDir,SessionOpt,MiceOpt)


%% Automated 
% Load all data
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%Predefine
AllTimeLines = cell(length(MiceOpt),size(SessionOpt,2));
AllInputs = [];
% Within folders, look for 'VR'
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    SessionOpttmp = SessionOpt(midx,:);
    SessionOpttmp(cellfun(@isempty,SessionOpttmp)) = [];
    for sesidx=1:length(SessionOpttmp)
        thisdate = Dates4Mouse{find(~(cellfun(@isempty,cellfun(@(X) strfind(SessionOpttmp{sesidx},X),Dates4Mouse,'UniformOutput',0))))};
        thisses = strsplit(SessionOpttmp{sesidx},[thisdate '\']);
        thisses = thisses{end};
     
        TimelineFiles = dir(fullfile(DataDir,MiceOpt{midx},thisdate,thisses,'*Timeline.mat'));
        
        if isempty(TimelineFiles)
            continue %No Timeline, continue
        else
            flag = 1;
        end
        % Latest trial has all information you need
        for tridx = 1
            tmp = load(fullfile(TimelineFiles(tridx).folder,TimelineFiles(tridx).name));
            tmptl = tmp.Timeline.rawDAQData;
            tmpinputs = {tmp.Timeline.hw.inputs(:).name};
            tmpinputs = {'timestamps',tmpinputs{:}};
            if isempty(AllInputs)
                AllInputs = tmpinputs;
            else
                if any(~ismember(AllInputs,tmpinputs))
                    warning('Different inputs?!')
                    keyboard
                end
            end
            tmpdat = cat(2,tmp.Timeline.rawDAQTimestamps',tmp.Timeline.rawDAQData);
            AllTimeLines{midx,sesidx} =  tmpdat;
            
        end
        
    end
    
end
