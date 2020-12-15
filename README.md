# CVMHT
Dataset and Code of CVMHT (Complementary-View Multiple Human Tracking), accepted in AAAI 2020.

```
@inproceedings{han2020cvmht,
  title={Complementary-View Multiple Human Tracking}, 
  author={Han, Ruize and Feng, Wei and Zhao, Jiewen and Niu, Zicheng and Zhang, Yunjun and Wan, Liang and Wang, Song},  
  year={2020},  
  booktitle={AAAI Conference on Artificial Intelligence}
}
```
## Introduction

<div align=center><img src="https://github.com/HanRuize/CVMHT/blob/master/figs/example2.jpg" width="500" height="150" alt="example"/><br/>
<div align= justify>

While tracking has many applications, typically two pieces of information can be provided by tracking results: accurate trajectories and appearance of the targets over time. This clearly introduces a conflict – if the camera is too close to the targets, limited coverage and frequent mutual occlusions prevent the accurate detection of their trajectories; if the camera is too far away from the targets, it is difficult to capture the detailed appearance of targets that are important for many applications such as person identification, action recognition, etc. 

In this paper, we present a new camera setting to address this problem. To track a group of people, which we refer to as subjects in this paper, on the ground, we use two cameras with different views and synchronized clock: A top-view camera at a high altitude, e.g, mounted to flying drone, provides a global birds-eye view of the subjects and the whole scene as shown in above Figure (a). A horizontal-view camera on the ground, e.g., mounted to a helmet worn by one person, which is static or moves/rotates smoothly without drastic visual field changes, captures the detailed appearance of subjects of interest, as shown in above Figure (b).

## Method

<div align=center><img src="https://github.com/HanRuize/CVMHT/blob/master/figs/solution2.jpg" width="500" height="250" alt="example"/><br/>
<div align= justify>
Given a pair of temporally-aligned videos that are taken from the top view and horizontal view, respectively, we first synchronously split these two videos into short clips with the same length. In each video clip, we extract a set of subject tracklets using a simple overlap criteria: detected subjects (in the form of bounding boxes) with good overlap between two adjacent frames, are connected to form tracklets. We refer to the resulting tracklets as single-view tracklets since they are extracted from the top-view and horizontal-view videos, respectively. We then conduct cross-view cross-clip data association for these single-view tracklets. More specifically, as shown in above Figure, the single-view tracklets from two adjacent clips in two views are fed into corresponding cluster as nodes. We then establish the subject association between clips and across views by using a joint optimization function. We refer to the generated tracking trajectories between two clips and across two views as cross-view short trajectories. Finally, we stitch the short trajectories by considering the frame overlap over time and obtain the cross-view long trajectories as the final tracking results.

## Experiments

![res](https://github.com/HanRuize/CVMHT/blob/master/figs/fig_case.png)

In horizontal-view videos, it is common to have subjects with mutual occlusion and being out-of-view. In this case, existing online trackers, e.g., DMAN can not associate the long-term lost subjects when they reappear in the view. Two examples are shown in above Figure. The top two rows show the case of mutual occlusions. From the top view at frame #180, we can find that two subjects (ID number 2, 3) are occluded by others and DMAN switches the ID of them when they reappear in the filed of view at frame #210. Our method keeps the original ID number. Similarly, we focus on the key subject (ID number 4) which goes out of view at frame #165 in the horizontal view. We can find that this subject is reassigned to a new ID number by DMAN. Our approach gets the original ID number of the target, which is consistent to its ID number in the top view.


Dataset: Link: https://pan.baidu.com/s/1VW9Ubr6AnBtn-F7U6THNzA Password: xnj9.

To get the annotation, please contact the authors. The dataset is only used for academic research.

Code: Mainly by Ruize Han (han_ruize@tju.edu.cn); Jiewen Zhao (zhaojw@tju.edu.cn).
