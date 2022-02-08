savefolder = fullfile('E:\Data\Histology\Sequence');
opts = dir(fullfile('E:\Data\Histology\CB007'))
if ~isdir(savefolder)
    mkdir(savefolder)
end
for i = 1:length(opts)
    
        imgs = dir(fullfile(opts(i).folder,opts(i).name,'*.tif'));
        if ~isempty(imgs)
            cellfun(@(X) copyfile(fullfile(X.folder,X.name),fullfile(savefolder,X.name)),{imgs},'UniformOutput',0)
        end
end