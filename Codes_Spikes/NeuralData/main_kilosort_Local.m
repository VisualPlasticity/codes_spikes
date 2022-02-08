%% you need to change most of the paths in this block
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab')) % for converting to Phy
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\Kilosort2.5'))
pathToYourConfigFile = 'C:\Users\EnnyB\Documents\MATLAB\LocalCode\'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'neuropixPhase3B2_kilosortChanMap.mat'; %3b2 used by Anwar, took same

%% Automated
% Load all data
% Fi nd available datasets (always using dates as folders)
DateOpt = cellfun(@(X) dir(fullfile(DataDir,X,'*-*')),MiceOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

for midx = 1 :length(MiceOpt)
    Dates4Mouse = DateOpt{midx};
    SessionOpttmp = SessionOpt(midx,:);
    SessionOpttmp(cellfun(@isempty,SessionOpttmp)) = [];
    for sesidx = 1:length(SessionOpttmp)
        
        thisdate = Dates4Mouse{find(~(cellfun(@isempty,cellfun(@(X) strfind(SessionOpttmp{sesidx},X),Dates4Mouse,'UniformOutput',0))))};
        thisses = strsplit(SessionOpttmp{sesidx},[thisdate '\']);
        thisses =thisses{end};

        rootS = fullfile(LocalDir,MiceOpt{midx},thisdate);
        if exist(rootS,'dir')
            disp([rootS 'already preprocessed by KiloSort, continue...'])
            continue
        else
            disp([rootS 'processing in KiloSort...'])
        end

        
        rootZ = fullfile('\\znas\Subjects\',MiceOpt{midx},thisdate,'ephys');
        tmpfoldx = dir(fullfile(rootZ,'*','*imec0'));
      

        if isempty(tmpfoldx)
            disp(['No Ephys data in ' rootZ ', check?'])
            continue
        end
        if any(~cellfun(@isempty,cellfun(@(X) strfind(X,'Concat'),{tmpfoldx(:).name},'UniformOutput',0)))
            idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,'Concat'),{tmpfoldx(:).name},'UniformOutput',0)));
            if exist(fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0.ap.bin']))
                tmpfoldx = tmpfoldx(idx);
                disp('Take concatenated data as ephys data')
            else
                tmpfoldx(idx)=[];
            end
        end
        if length(tmpfoldx)>1
            disp('Found multiple ephys recordings. Add together?')
            flag = 0;
            while ~flag
                S = input('(Y/N)','s');
                if strcmpi(S,'Y')
                    disp('Adding together ... ')
                    mkdir(fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat']))
                    mkdir(fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0']))
                    fileout = fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0.ap.bin'])
                    
                    fout = fopen(fileout,'w');
                    for idx = 1:length(tmpfoldx)
                        % find the binary file
                        fs          = [dir(fullfile(tmpfoldx(idx).folder,tmpfoldx(idx).name, '*ap.bin')) dir(fullfile(tmpfoldx(idx).folder,tmpfoldx(idx).name, '*.dat'))];
                        fin = fopen(fullfile(fs(1).folder,fs(1).name),'r');
                        bytes       = get_file_size(fullfile(fs(1).folder,fs(1).name)); % size in bytes of raw binary
                        
                        while ~feof(fin)
                            temp = fread(fin,10000,'int16');
                            fwrite(fout,temp,'int16');
                            fprintf(fout,'\n');
                        end
                        fclose(fin)
                    end
                    fclose(fout)
                    
                    % Also add the nidq file together
                    fileout = fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.nidq.bin'])
                    
                    fout = fopen(fileout,'w');
                    for idx = 1:length(tmpfoldx)
                        % find the binary file
                        fs          = dir(fullfile(tmpfoldx(idx).folder, '*.bin'));
                        fin = fopen(fullfile(fs(1).folder,fs(1).name),'r');
                        bytes       = get_file_size(fullfile(fs(1).folder,fs(1).name)); % size in bytes of raw binary
                        
                        while ~feof(fin)
                            temp = fread(fin,10000,'int16');
                            fwrite(fout,temp,'int16');
                            fprintf(fout,'\n');
                        end
                        fclose(fin)
                    end
                    fclose(fout)
                    
                    % and the lfp
                     fileout = fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0'],[thisdate '_' MiceOpt{midx} '_g0_Concat_t0.imec0.lf.bin'])
                    
                    fout = fopen(fileout,'w');
                    for idx = 1:length(tmpfoldx)
                        % find the binary file
                        fs          = [dir(fullfile(tmpfoldx(idx).folder,tmpfoldx(idx).name, '*lf.bin')) dir(fullfile(tmpfoldx(idx).folder,tmpfoldx(idx).name, '*.dat'))];
                        fin = fopen(fullfile(fs(1).folder,fs(1).name),'r');
                        bytes       = get_file_size(fullfile(fs(1).folder,fs(1).name)); % size in bytes of raw binary
                        
                        while ~feof(fin)
                            temp = fread(fin,10000,'int16');
                            fwrite(fout,temp,'int16');
                            fprintf(fout,'\n');
                        end
                        fclose(fin)
                    end
                    fclose(fout)
                    
                    tmpfoldx = dir(fullfile(rootZ,[thisdate '_' MiceOpt{midx} '_g0_Concat'],'*imec0'));
                    flag=1;
                elseif strcmp(S,'N')
                    disp('Not adding together ... ')
                    disp('Which one do you want to have processed then (choose main directory)?')
                    tmp = uigetdir(rootZ);
                    tmpfoldx = dir(fullfile(tmp,'*imec0'));
                    flag = 1;
                else
                    disp('Follow instructions please')
                end
            end
        end
    
            mkdir(rootS)

        rootZ = fullfile(tmpfoldx(1).folder,tmpfoldx(1).name); % the raw data binary file is in this folder
        rootH = fullfile(tmpdatafolder,MiceOpt{midx},thisdate);% path to temporary binary file (same size as data, should be on fast SSD)
        ops.trange    = [0 Inf]; % time range to sort
        ops.NchanTOT  = 385; % total number of channels in your recording
        
        run(fullfile(pathToYourConfigFile, 'configFile384.m'))
        ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
        ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);
        
        %% this block runs all the steps of the algorithm
        fprintf('Looking for data inside %s \n', rootZ)
        
        % make temporary directory if nonexistent
        if ~isdir(rootH)
            mkdir(rootH)
        end
        % main parameter changes from Kilosort2 to v2.5
        ops.sig        = 20;  % spatial smoothness constant for registration
        ops.fshigh     = 300; % high-pass more aggresively
        ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
        
        % is there a channel map file in this folder?
        fs = dir(fullfile(rootZ, 'chan*.mat'));
        if ~isempty(fs)
            ops.chanMap = fullfile(rootZ, fs(1).name);
        end
        
        % find the binary file
        fs          = [dir(fullfile(rootZ, '*.bin'))' dir(fullfile(rootZ, '*.dat'))];
        ops.fbinary = fullfile(rootZ, fs(1).name);
        
        % preprocess data to create temp_wh.dat
        rez = preprocessDataSub(ops);
        %
        % NEW STEP TO DO DATA REGISTRATION
        rez = datashift2(rez, 1); % last input is for shifting data
        
        % ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
        iseed = 1;
        
        % main tracking and template matching algorithm
        rez = learnAndSolve8b(rez, iseed);
        
        % OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
        % See issue 29: https://github.com/MouseLand/Kilosort/issues/29
        rez = remove_ks2_duplicate_spikes(rez);
        
        % final merges
        rez = find_merges(rez, 1);
        
        % final splits by SVD
        rez = splitAllClusters(rez, 1);
        
        % decide on cutoff
        rez = set_cutoff(rez);
        % eliminate widely spread waveforms (likely noise)
        rez.good = get_good_units(rez);
        
        fprintf('found %d good units \n', sum(rez.good>0))
        
        % write to Phy
        fprintf('Saving results to Phy  \n')
        rezToPhy(rez, rootS);
        
        %% if you want to save the results to a Matlab file...
        if 0
            % discard features in final rez file (too slow to save)
            rez.cProj = [];
            rez.cProjPC = [];
            
            % final time sorting of spikes, for apps that use st3 directly
            [~, isort]   = sortrows(rez.st3);
            rez.st3      = rez.st3(isort, :);
            
            % Ensure all GPU arrays are transferred to CPU side before saving to .mat
            rez_fields = fieldnames(rez);
            for i = 1:numel(rez_fields)
                field_name = rez_fields{i};
                if(isa(rez.(field_name), 'gpuArray'))
                    rez.(field_name) = gather(rez.(field_name));
                end
            end
            
            % save final results as rez2
            fprintf('Saving final results in rez2  \n')
            fname = fullfile(rootZ, 'rez2.mat');
            save(fname, 'rez', '-v7.3');
        end
    end
end