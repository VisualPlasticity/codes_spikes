% calculate Contrast for VR screenshots

% which files?
scrshtdir = 'X:\Enny\fUSI\ScreenShotsCorridor'
scrshtlist = dir(fullfile(scrshtdir,'*.jpg'));

% find options for contrast and position
contrast=cellfun(@(X) str2num(X{1}),cellfun(@(X) strsplit(X{2},'_'),cellfun(@(X) strsplit(X,'Contrast_'),{scrshtlist(:).name},'UniformOutput',0),'UniformOutput',0),'UniformOutput',0);
ContrOpt = unique(cell2mat(contrast));
position=cellfun(@(X) str2num(X{1}),cellfun(@(X) strsplit(X{2},'.jpg'),cellfun(@(X) strsplit(X,'Pos_'),{scrshtlist(:).name},'UniformOutput',0),'UniformOutput',0),'UniformOutput',0);
PosOpt = unique(cell2mat(position));

% Prepare RMS
RMS = nan(length(PosOpt),length(ContrOpt));

%Define contrast
for i = 1:length(scrshtlist)
    disp([num2str(i) ' out of ' num2str(length(scrshtlist))])
    %What contrast, what position?
    nameparts=strsplit(scrshtlist(i).name,'Contrast_');
    nameparts = strsplit(nameparts{2},'_');
    thiscontrid = find(ismember(ContrOpt,str2num(nameparts{1})));
    nameparts = strsplit(nameparts{3},'.jpg');
    thispos = find(ismember(PosOpt,str2num(nameparts{1})));
    
    %Read image 
    Img=imread(fullfile(scrshtlist(i).folder,scrshtlist(i).name));
    Img = double(rgb2gray(Img))./255; %make it intensity
    
    %Remove photodiode square
    Img(1:190,1:170)=nan; %Check position
    if i==1
        figure; imagesc(Img)
        drawnow
    end
     
    % Average
    AvgIntensity = nanmean(Img(:)); %faster to calculate this once then doing this in the below function
    %Root Mean Square contrast
    PixelSq = arrayfun(@(X) (X-AvgIntensity).^2,Img);  
    
    % RMS formula
    RMS(thispos,thiscontrid) = sqrt(nansum(PixelSq(:))/numel(Img));    
end

figure; 
imagesc(RMS')
colormap gray
set(gca,'YTick',1:length(ContrOpt),'YTickLabel',arrayfun(@(X) num2str(X),ContrOpt,'UniformOutput',0),'XTick',1:5:length(PosOpt),...
    'XTickLabel',arrayfun(@(X) num2str(X),PosOpt(1:5:end),'UniformOutput',0),'YDir','normal')
ylabel('ContrastSetting')
xlabel('Position')
 
box off
colorbar
title('RMS Contrast')

saveas(gcf,'C:\Users\EnnyB\Dropbox\UCL_London\SetUp\ContrastAcrossCorridor.fig')
saveas(gcf,'C:\Users\EnnyB\Dropbox\UCL_London\SetUp\ContrastAcrossCorridor.bmp')
save('C:\Users\EnnyB\Dropbox\UCL_London\SetUp\ContrastAcrossCorridor.mat','RMS','PosOpt','ContrOpt')
