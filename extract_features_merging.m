function [tv_midLevelTracklet,hv_midLevelTracklet] = extract_features_merging(NN,  NN_original, nodes, cv_segments, seqName,param_tracklet, flag_visualize_midLevelTracklets)

saveFeat = ['./Data/Features/' seqName '_merging_feat/'];
if(~exist(saveFeat))
    mkdir(saveFeat);
end
%% Sequence Info
data_directory = './Data/Images/';
imPath = fullfile(data_directory,seqName,'/');
images = dir([imPath '/*.jpg']);
numOfClusters = param_tracklet.num_cluster; 
segmentLength = param_tracklet.num_frames; 

for iFile = 1:length(NN)   
    
    NN_curr = NN{iFile};
    nodes_curr = nodes{iFile};
    NN_original_curr = NN_original{iFile};
    
    segment_num = (iFile-1)*round((numOfClusters)/2);
    numClicks = size(nodes_curr,1);
    cc = hsv(numClicks);
    cumulativeSUM = cumsum(NN_curr);
    centers = cell(1,numClicks);

    if (flag_visualize_midLevelTracklets)
        im1 = imread([imPath images((segment_num)*segmentLength+1).name]);
        im2 = imread([imPath images(min((segment_num+numOfClusters)*segmentLength+1,length(images))).name]);
        im = (im2double(im1)+im2double(im2))/2;
        imshow(im);
        hold on;
    end
    
    tv_tracklet_cnt=1;
    hv_tracklet_cnt=1;
    dummy_tracklet_flag = zeros(1,numClicks);
    tv_tracklet_flag = zeros(1,numClicks);
    hv_tracklet_flag = zeros(1,numClicks);
    for iClick = 1:numClicks     
        tv_detections = [];
        tv_frames = [];
        tv_tracklet_color_hist = [];
        hv_detections = [];
        hv_frames = [];
        hv_tracklet_color_hist = [];
        cnt_dummy_nodes = 0;
        for ii = [1,2,4,3] %1:size(nodes_curr,2)
            if(ii==1)
                trackletID = nodes_curr(iClick,ii);
            else
                trackletID = nodes_curr(iClick,ii)-cumulativeSUM(ii-1);
            end
            if(trackletID>NN_original_curr(ii))
                cnt_dummy_nodes = cnt_dummy_nodes+1;
                continue;
            end
            seg_i = [1,2,2,1];                
            iii = seg_i(ii);
            if ii <= 2  % Top-view middle-level tracklet
                segment = cv_segments{1}; 
                try
                tv_detections = [tv_detections ; segment{segment_num+iii}.tracklet(trackletID).detection];
                tv_frames = [tv_frames ; segment{segment_num+iii}.tracklet(trackletID).frame];
                catch err
                end
                
                aa = segment{segment_num+iii}.tracklet.color_hist_median;
                tv_tracklet_color_hist = [tv_tracklet_color_hist aa(:)];
            else  % Hor-view middle-level tracklet
                segment = cv_segments{2};
                hv_detections = [hv_detections ; segment{segment_num+iii}.tracklet(trackletID).detection];                
                hv_frames = [hv_frames ; segment{segment_num+iii}.tracklet(trackletID).frame];
                aa = segment{segment_num+iii}.tracklet.color_hist_median;
                hv_tracklet_color_hist = [hv_tracklet_color_hist aa(:)];              
            end         
        end
        
        if (isempty(tv_detections)&&isempty(hv_detections))
            continue;
        end
        
        % We remove nodes which have only one real node
        if(cnt_dummy_nodes< (size(nodes_curr,2)-1))
            % Store the tracklet information
            if ~isempty(tv_detections)
                tv_tracklet.color_hist{tv_tracklet_cnt} = tv_tracklet_color_hist;
                tv_tracklet.bbox{tv_tracklet_cnt}       = tv_detections;
                tv_tracklet.frames{tv_tracklet_cnt}     = tv_frames; 
                tv_tracklet_flag(1,iClick) = 1;
                tv_tracklet_cnt = tv_tracklet_cnt+1;              
            end
            if ~isempty(hv_detections)          
                hv_tracklet.color_hist{hv_tracklet_cnt} = hv_tracklet_color_hist;
                hv_tracklet.bbox{hv_tracklet_cnt}       = hv_detections;
                hv_tracklet.frames{hv_tracklet_cnt}     = hv_frames;
                hv_tracklet_flag(1,iClick) = 1;
                hv_tracklet_cnt = hv_tracklet_cnt+1;
            end
                        
            if (flag_visualize_midLevelTracklets)                
                detections = hv_detections;
                try
                centers{iClick}=[(detections(:,1)+detections(:,3))/2 (detections(:,2)+detections(:,4))/2];
                catch err
                end
                plot(centers{iClick}(:,1),centers{iClick}(:,2),'-mo','color',cc(iClick,:),'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',cc(iClick,:));
                pause(0.4);
                hold on;
            end
        end
    end

    tv_midLevelTracklet(iFile).flag =  tv_tracklet_flag;
    hv_midLevelTracklet(iFile).flag =  hv_tracklet_flag;
    try
    tv_midLevelTracklet(iFile).color_hist = tv_tracklet.color_hist;
    catch err
    end
    tv_midLevelTracklet(iFile).bbox       = tv_tracklet.bbox;
    tv_midLevelTracklet(iFile).frames     = tv_tracklet.frames;
    hv_midLevelTracklet(iFile).color_hist = hv_tracklet.color_hist;
    hv_midLevelTracklet(iFile).bbox       = hv_tracklet.bbox;
    hv_midLevelTracklet(iFile).frames     = hv_tracklet.frames;
    for i = 1 : length(hv_tracklet.frames)
    if min(hv_tracklet.frames{i}) < (iFile - 1)*10
         hv_tracklet.frames;           
    end
    end

    %save([saveFeat sprintf('segment_%03d.mat',iFile)],'tracklet');
    clear tracklet hv_tracklet tv_tracklet;
end

