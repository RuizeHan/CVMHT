%% CVMOT
%% 07/2019
%% Ruize Han, Jiewen Zhao

% Based on
%% GMMCP Tracker
%% Afshin Dehghan, Shayan Modiri Assari

clc;clear;close all;
warning off;
addpath('./MOT_toolbox/');
addpath('./toolbox/')
addpath('./MHA/')
addpath('./MHA/MHA_code')
addpath('./MHA/MHA_code/unit')
% addpath('./MOT_toolbox/cplex/');
% addpath('./MOT_toolbox/cplexWindows/')
addpath('./MOT_toolbox/GOG/')
addpath('/opt/ibm/ILOG/CPLEX_Studio129/cplex/matlab/x86-64_linux/');  % todo : Change to the installation path of CPLEX_Studio

cd '/media/han/VPG_T4_#1/Raiser/AAAI20/CVMOT_Tracker_github'  % todo : Change to the path of the CVMHT code 
%% Sequence Info
seq = configSeqs;

% Parameter setting
sv_appearanceWeight = 0.5;
Aweight = 0.3; 
dummyWeight = 0.3;
sv_netCost_norm = 1;
cv_netCost_norm = 0;
occ_angle = 2;
motion_punish = 0; 
spatial_punish = 0;

param.motion_punish = motion_punish;
param_cv_netCost.spatial_punish = spatial_punish;
params.dummyWeight = dummyWeight;
    
Res_name = 'Res_tracking';

tracklet_time_seq = zeros(length(seq),1);
sv_ass_time_seq = zeros(length(seq),1);
cv_ass_time_seq = zeros(length(seq),1);
solution_time_seq = zeros(length(seq),1);

for seq_i = 1  : length(seq)
scene_name =  seq{seq_i}.name; 

view_names = {'top','hor'};
params.det_type = 2;    % 1:detection 2:GT
params.search_gopro = 2; % 1: auto, 0:GT  2:Tracking res;
params.load_detection = 1;
params.load_svcost = 1;
params.load_cvcost = 1;
params.load_solusion = 1;
params.bx_min_size = [0,3000];
params.vis_tracklet = 0;
params.vis_result = 1;
params.sv_netCost_norm = sv_netCost_norm;
params.cv_netCost_norm = cv_netCost_norm;

params.sv_appearanceWeight = sv_appearanceWeight;
param_cv_netCost.Aweight = Aweight;
param_cv_netCost.Sweight = 1 - param_cv_netCost.Aweight;
param.overlap_rho = 1.0;

flag_visualization_ll_tracklet = 0;
flag_visualize_midLevelTracklets = 0;

data.cv_segments = cell(2,1);
data.cv_detections = cell(2,1);
data.cv_images = cell(2,1);
data.cv_im_directory = cell(2,1);

for view_i =  1 : 2
view = view_names{view_i};
sequence_name = [scene_name,'_',view];
dataset_directory = '../dataset';      % todo : Change to the path of the CVMHT dataset 
data_directory = './Data/';
im_directory = fullfile(dataset_directory,'Images/',scene_name,'frame_sel/',view);
images = dir(fullfile(im_directory,'*.jpg'));

data.cv_im_directory{view_i} = im_directory;
data.cv_images{view_i} = images;

% % load the annotion bounding boxes
load(fullfile(dataset_directory,'GT',[sequence_name,'.mat']));
GT = objs_inter(:,2)';
    
if params.det_type == 2
    detections = GT;
else 
    % load the detection bounding boxes
    load(fullfile(dataset_directory,'Detection',[sequence_name,'.mat']));
    detections = det_objs(:,2)';
    % detections(1:10) = GT(1:10);
end


% ----- visualize the bounding box -------
% image = images(1);
% img = imread([im_directory,'/',image.name]);
% det = detections{1};
% imshow(img);
% for i = 1 : size(det,1)
% det1 = det(i,:);
% rec1(:,1) = det1(:,1);
% rec1(:,2) = det1(:,2);
% rec1(:,3) = det1(:,3) - det1(:,1);
% rec1(:,4) = det1(:,4) - det1(:,2);
% 
% rectangle('Position',rec1(1,1:4),'LineWidth',2,'LineStyle','-','EdgeColor','r');hold on;
% end


%----- detection transformation for Top&Hor Dataset --------
detections_bx_filter = cell(length(detections),1);
detections_bx = cell(length(detections),1);
detections_filter = cell(length(detections),1);
for i = 1:length(detections)
    
  if params.det_type == 2 % the GT labels
      detections_bx{i}(:,5) = detections{i}(:,5);  
  else
      detections{i}(:,5) = 1;
  end
  
  % filter the detection with too small size
  detections_bx_w = detections{i}(:,3) - detections{i}(:,1);
  detections_bx_h = detections{i}(:,4) - detections{i}(:,2);
  filter_index = detections_bx_w .* detections_bx_h > ones(length(detections_bx_w),1) * params.bx_min_size(view_i);
  
  detections_bx_filter{i}(:,1) = detections{i}(filter_index,1);
  detections_bx_filter{i}(:,2) = detections{i}(filter_index,2);
  detections_bx_filter{i}(:,3) = detections{i}(filter_index,3) - detections{i}(filter_index,1);
  detections_bx_filter{i}(:,4) = detections{i}(filter_index,4) - detections{i}(filter_index,2); 
  detections_bx_filter{i}(:,5) = detections{i}(filter_index,5); 
  detections_filter{i}(:,1:5) = detections{i}(filter_index,1:5);

end

data.cv_detections{view_i} = detections_bx_filter;

%% Initialize the parameters
[param_tracklet,param_merging,param_tracking,param_netCost]=set_param_gmmcp;
%  MHA parameters
param_tracklet.Lambda = 0.015;
param_tracklet.Rho = 15;
param_tracklet.occ_angle = occ_angle;

param_tracklet.seqName = sequence_name;
param_tracklet.sceneName = scene_name;
param_tracklet.search_gopro = params.search_gopro;
param_tracklet.rho = param.overlap_rho;
param_tracklet.data_directory = data_directory;
param_tracklet.dataset_directory = dataset_directory;
param_tracklet.num_segment = round(length(images)/param_tracklet.num_frames);

param_netCost.motion_punish = param.motion_punish;
param_netCost.seqName = sequence_name;


%% Create Low-Level Tracklets and Extract Appearance Features
if(exist(fullfile(param_tracklet.data_directory,'Features',['tracklets_' param_tracklet.seqName '_nf.mat']), 'file') && params.load_detection == 1 )
    load(fullfile(param_tracklet.data_directory,'Features',['tracklets_' param_tracklet.seqName '_nf.mat']));
else
    fprintf('Creating Low-level Tracklets  / ');
    tt = tic;
    segment = ll_tracklet_generator_nf(im_directory,detections_filter,param_tracklet, images, flag_visualization_ll_tracklet);
    tracklet_time_seq(seq_i) = tracklet_time_seq(seq_i) + toc(tt);
    fprintf('\nTime Elapsed:%0.2f\n',toc(tt));
end
data.cv_segments{view_i} = segment;
end

param_tracklet.num_segment = length(segment);

%% Create NetCost Matrix
fprintf('Creating SV-NetCost Matrix for Low-level Tracklets  / ');
tt = tic;

if(exist(fullfile(param_tracklet.data_directory,'Cost',['sv_net_cost_' param_tracklet.seqName '.mat']), 'file') && params.load_svcost == 1)
    load(fullfile(param_tracklet.data_directory,'Cost',['sv_net_cost_' param_tracklet.seqName '.mat']));
else
    sv_net_cost = cell(2,1); 
    for view_i = 1 : 2
        net_cost = create_netCost(data.cv_segments{view_i},param_netCost);
        sv_net_cost{view_i} = net_cost;
    end
    sv_ass_time_seq(seq_i) = sv_ass_time_seq(seq_i) + toc(tt);
    save(fullfile(param_tracklet.data_directory,'Cost',['sv_net_cost_' param_tracklet.seqName '.mat']),'sv_net_cost','-v7.3');
end

fprintf('Creating CV-NetCost Matrix for Low-level Tracklets  / ');
tt = tic;
if(exist(fullfile(param_tracklet.data_directory,'Cost',['cv_net_cost_' param_tracklet.seqName '.mat']), 'file') && params.load_cvcost == 1)
    load(fullfile(param_tracklet.data_directory,'Cost',['cv_net_cost_' param_tracklet.seqName '.mat']));
else
    cv_net_cost = creat_cv_netCost(data,param_cv_netCost,param_tracklet);  
    cv_ass_time_seq(seq_i) = cv_ass_time_seq(seq_i) + toc(tt);
    save(fullfile(param_tracklet.data_directory,'Cost',['cv_net_cost_' param_tracklet.seqName '.mat']),'cv_net_cost','-v7.3');
end

%% 
% ----- visualize the tracklet bounding box -------
if params.vis_tracklet == 1
for seg_k = 31 % 1 : param_tracklet.num_segment
    for seg_i = seg_k : seg_k + 1
    for view_i =  1 : 2       
        im_directory = data.cv_im_directory{view_i};
        images = data.cv_images{view_i};         
        segments = data.cv_segments{view_i};
        segment = segments{seg_i};
        iImg = seg_i * 10;
        image = images(iImg);
        img = imread([im_directory,'/',image.name]);
         figure
         imshow(img); hold on; 
         set (gcf,'Position',[600 * (view_i-1) ,500 * (seg_i-seg_k),600,400])
         tracklet = segment.tracklet;
         for i=1:size(tracklet,2)    
             % c = hsv(size(tracks,2));    
             det = tracklet(i).detection;
             xc = round((det(:,1)+ det(:,3))/2);
             yc = round((det(:,2)+ det(:,4))/2);
             % plot(xc,yc,'-o','color','r','MarkerFaceColor','r','LineWidth',1);
             text(xc(1),yc(1),num2str(i),'color','g','fontsize',15);      
         end     
    end
    end
%     close all
end
end

fprintf('\nTime Elapsed:%0.2f\n',toc(tt));
%% 
tt = tic;
%% Run GMMCP with ADN on non-overlaping segments
filename = {'nodes','NN','NN_original'};
if(exist(fullfile(param_tracklet.data_directory,'GMMCP_Tracklet',['nodes_' param_tracklet.seqName '.mat']), 'file') && params.load_solusion == 1)
    for f = 1:length(filename)
        load(fullfile(param_tracklet.data_directory,'GMMCP_Tracklet',[filename{f} '_' param_tracklet.seqName '.mat']));
    end
else
    cnt_batch=1;
    NN = cell(1,param_tracklet.num_segment-1);
    NN_original = cell(1,param_tracklet.num_segment-1);
    nodes = cell(1,param_tracklet.num_segment-1);
    for iSegment=1:round(param_tracklet.num_cluster/2):param_tracklet.num_segment
        if iSegment < param_tracklet.num_segment
           fprintf('computing tracklets for segment %d to %d \n',iSegment,min(iSegment + param_tracklet.num_cluster-1,param_tracklet.num_segment));
           [NN{cnt_batch},NN_original{cnt_batch}, nodes{cnt_batch}] = GMMCP_Tracklet_Generation(sv_net_cost,cv_net_cost, iSegment,...
             min(iSegment+param_tracklet.num_cluster-1,param_tracklet.num_segment),params);
            cnt_batch = cnt_batch+1;
        end
    end
    for f = 1:length(filename)
      save(fullfile(param_tracklet.data_directory,'GMMCP_Tracklet',[filename{f} '_' param_tracklet.seqName '.mat']),filename{f},'-v7.3');
    end
end

%% create NetCost Matrix for Merging
cv_midLevelTracklets = cell(2,1);
% sequence_name = view_names{2};
% images = data.cv_images{2};
[cv_midLevelTracklets{1},cv_midLevelTracklets{2}] = extract_features_merging(NN, NN_original, nodes, data.cv_segments, sequence_name,param_tracklet, flag_visualize_midLevelTracklets);


%% Stitch the tracklets (Final Data Association)
fprintf('Stitching Tracklets to form final tracks \n ');

finalTracks = cell(2,1);
trackRes = cell(2,1); 

[cv_midLevelTracklets, trackRes] = stitchTracklets(cv_midLevelTracklets);
solution_time_seq(seq_i) = solution_time_seq(seq_i) + toc(tt);

if(~exist([param_tracklet.data_directory,Res_name],'dir'))
    mkdir([param_tracklet.data_directory,Res_name]);
end

tv_res = trackRes{1};
hv_res = trackRes{2};

save(fullfile(param_tracklet.data_directory,Res_name,['track_res_' param_tracklet.sceneName '_top.mat']),'tv_res','-v7.3');
save(fullfile(param_tracklet.data_directory,Res_name,['track_res_' param_tracklet.sceneName '_hor.mat']),'hv_res','-v7.3');

fprintf('\nTime Elapsed:%0.2f\n',toc(tt));
% Visualize Final Trackig Results
if params.vis_result == 1
    outDir = cell(2,1);
    for view_i = 1 : 2  
        outDir{view_i} = sprintf('./Data/vis_result/%s/%s/%s/', Res_name, param_tracklet.sceneName, view_names{view_i});
%       outDir{view_i} = sprintf('./trackingResults/%s/',view_names{view_i});
    end
    plotTracking(trackRes, data.cv_im_directory,data.cv_images,0,outDir);
end

end

