function Listvideofiles(DataDir,MiceOpt,storevideopath,copybodytolocal)

%% Automated
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

disp('List of videos for eyetracking to run in python deeplabcut:')
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        SessionOpttmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'*10*'));%All my sessions start with 10(1/2/3/etc.)        
        for sesidx=1:length(SessionOpttmp)
            % Al Video data
            EyeVideoFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*eye.mj2*'));
            if ~isempty(EyeVideoFiles) 
                namecomp = strsplit(EyeVideoFiles(:).name,'.mj2');
                csvfiles = dir(fullfile(EyeVideoFiles(:).folder,[namecomp{1} '*.csv']));
                if isempty(csvfiles)
                    disp(['r''', fullfile(EyeVideoFiles(:).folder,EyeVideoFiles(:).name) ','])                   
                end
            end
        end
    end
end

disp('List of videos for belly to run in python deeplabcut:')
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        SessionOpttmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'*10*'));%All my sessions start with 10(1/2/3/etc.)
        
        for sesidx=1:length(SessionOpttmp)
            BellyVideoFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*bellyCam.mj2'));
            if ~exist(fullfile(storevideopath,'belly',BellyVideoFiles(:).name)) && ~isempty(BellyVideoFiles)
                disp(['r''', fullfile(BellyVideoFiles(:).folder,BellyVideoFiles(:).name) ','])
                if copybodytolocal
                    if ~isdir(fullfile(storevideopath,'belly'))
                        mkdir(fullfile(storevideopath,'belly'))
                    end
                    copyfile(fullfile(BellyVideoFiles(:).folder,BellyVideoFiles(:).name),fullfile(storevideopath,'belly',BellyVideoFiles(:).name));
                end
            end
          
        end
    end
end

disp('List of videos for body to run in python deeplabcut:')
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        SessionOpttmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'*10*'));%All my sessions start with 10(1/2/3/etc.)        
        for sesidx=1:length(SessionOpttmp)                     
            BodyideoFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*bodyCam.mj2s'));
            if ~exist(fullfile(storevideopath,'body',BodyideoFiles(:).name)) && ~isempty(BodyideoFiles)
                disp(['r''', fullfile(BodyideoFiles(:).folder,BodyideoFiles(:).name) ','])
                if copybodytolocal
                    if ~isdir(fullfile(storevideopath,'body'))
                        mkdir(fullfile(storevideopath,'body'))
                    end
                    copyfile(fullfile(BodyideoFiles(:).folder,BodyideoFiles(:).name),fullfile(storevideopath,'body',BodyideoFiles(:).name));
                end
            end
        end
    end
end
