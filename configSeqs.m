function frames = configSeqs


frame_sel = {
    struct('name','V0-S_canteen_0-G_1','startFrame',2,'endFrame',601,'topStart',2),... %602/601
    struct('name','V0-S_canteen_2-G_1','startFrame',2,'endFrame',601,'topStart',2),... %601/601
    struct('name','V0-S_canteen_0-G_3','startFrame',2,'endFrame',601,'topStart',2),...  %601/601
    struct('name','V0-S_canteen_1-G_3','startFrame',2,'endFrame',601,'topStart',2),... %601/601
    struct('name','V0-S_canteen_2-G_3','startFrame',2,'endFrame',601,'topStart',2),...  %601/601
    struct('name','V1-S_55-G_1','startFrame',601,'endFrame',1800,'topStart',601),...  %2701/2701
    struct('name','V1-S_55-G_2','startFrame',601,'endFrame',1500,'topStart',601),...  %2701/2701
    struct('name','V1-S_library-G_1','startFrame',1,'endFrame',1200,'topStart',1),...  %2700/2700
    struct('name','V1-S_library-G_2','startFrame',1,'endFrame',900,'topStart',1),...  %2700/2700
    struct('name','V1-S_library-G_3','startFrame',1,'endFrame',600,'topStart',1),...  %2700/2700
    struct('name','V1-S_square-G_1','startFrame',1,'endFrame',1200,'topStart',1),...  %2701/2701
    struct('name','V1-S_square-G_2','startFrame',1,'endFrame',900,'topStart',1),...  %2701/2701
    struct('name','V1-S_square-G_3','startFrame',601,'endFrame',1200,'topStart',1),...  %2701/2701
    struct('name','V2-S_playground_1-G_2','startFrame',3601,'endFrame',4200,'topStart',3601),...  %4800/4800
    struct('name','V2-S_playground_2-G_3','startFrame',1801,'endFrame',2400,'topStart',1801),...  %4800/4800
};

%     struct('name','V0-S_55-G_1','startFrame',2,'endFrame',601),...  %602/601
%     struct('name','V0-S_55-G_2','startFrame',2,'endFrame',601),...  %602/601
%     struct('name','V0-S_55-G_2','startFrame',2,'endFrame',601),...  %602/601%   
%     struct('name','V1-S_55-G_3','startFrame',601,'endFrame',1200),...  %2701/2701  --out 
%     struct('name','V0-S_canteen_1-G_1','startFrame',2,'endFrame',601),... %601/627  --out  
frame_train = {
% trainging for top detection
%     struct('name','V1-S_55-G_1','startFrame',2501,'endFrame',2700),...  %2701/2701      
%     struct('name','V1-S_library-G_1','startFrame',2501,'endFrame',2700),...  %2700/2700
%     struct('name','V1-S_square-G_1','startFrame',2501,'endFrame',2700),...  %2701/2701   
% trainging for appearance feature
%     struct('name','V2-S_playground-G_2','startFrame',4501,'endFrame',4800),...  %4800/4800 
%     struct('name','V2-S_playground-G_3','startFrame',4501,'endFrame',4800),...  %4800/4800 
      struct('name','V0-S_canteen-G_3_5','startFrame',2,'endFrame',201),...  %601/601 % Manual annotation
      struct('name','V0-S_canteen-G_1_5','startFrame',2,'endFrame',201),...  %601/601 % Manual annotation
};
frames = frame_sel;
% frames = frame_train;
%     old
%     struct('name','V1-S_library-G_3','startFrame',1201,'endFrame',1800),...  %2700/2700
%     struct('name','V1-S_square-G_3','startFrame',1501,'endFrame',2100),...  %2701/2701