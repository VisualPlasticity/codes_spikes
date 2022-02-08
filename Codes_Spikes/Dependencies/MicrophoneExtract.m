function [AllMicrophone,AllFS,Allnbits] = MicrophoneExtract(DataDir,MiceOpt)


%% Automated 
% Load all data
% Find available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%Predefine
AllMicrophone = cell(length(MiceOpt),max(cell2mat(cellfun(@length,DateOpt,'UniformOutput',0))));
AllFS = AllMicrophone;
Allnbits = AllMicrophone;
% Within folders, look for 'VR'
for midx = 1:length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    sescount = 1;
    for didx = 1:length(Dates4Mouse)
        SessionOpttmp = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},'*10*'));%All my sessions start with 10(1/2/3/etc.)                
      
        for sesidx=1:length(SessionOpttmp)
            VRDataFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*VRBehavior*'));

            MicDataFiles = dir(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,'*mic*'));
            
            if isempty(VRDataFiles) || isempty(MicDataFiles)
                continue %Not VR and/or no recording found, continue
            else
                flag = 1;
            end
            % Latest trial has all information you need
            for tridx = 1
               tmp = load(fullfile(DataDir,MiceOpt{midx},Dates4Mouse{didx},SessionOpttmp(sesidx).name,MicDataFiles(tridx).name));
               
               if isfield(tmp,'micdata')
                   AllMicrophone{midx,sescount} = tmp.micdata;
               elseif isfield(tmp,'micData')
                   AllMicrophone{midx,sescount} = tmp.micData;
               else
                   keyboard
               end
              AllFS{midx,sescount} = tmp.Fs;
              Allnbits{midx,sescount} = tmp.nBits;
            end            
        end

            sescount = sescount+1;
    end

end
