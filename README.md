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

While tracking has many applications, typically two pieces of information can be provided by tracking results: accurate trajectories and appearance of the targets over time. This clearly introduces a conflict – if the camera is too close to the targets, limited coverage and frequent mutual occlusions prevent the accurate detection of their trajectories; if the camera is too far away from the targets, it is difficult to capture the detailed appearance of targets that are important for many applications such as person identification, ac-


An illustration of the top-view (a) and horizontalview (b) videos. The former is taken by a camera mounted to a drone in the air and the latter is taken by a GoPro worn by a wearer who walked on the ground. The proposed method jointly tracking multiple subjects, indicated by identical color boxes, across the two videos. Note that the global motion trajectory and local appearance are well presented in these
two complementary views.

## Method

![example](https://github.com/HanRuize/DSAR-CF/blob/master/figs/example.png)

Given a pair of temporally-aligned videos that are taken from the top view and horizontal view, respectively, we first synchronously split these two videos into short clips with the same length, e.g., 10 frames. In each video clip, we extract a set of subject tracklets using a simple overlap criteria: detected subjects (in the form of bounding boxes) with good overlap, e.g., higher than 50%, between two adjacent frames, are connected to form tracklets. We discard the tracklets that are overly short, e.g., traversing less than 3 frames. We refer to the resulting tracklets as single-view tracklets since they are extracted from the top-view and horizontal-view videos, respectively. We then conduct cross-view cross-clip data association for these single-view tracklets. More specifically, as shown in Fig. 2, the single-view tracklets from two adjacent clips in two views are fed into corresponding cluster

## Experiments

![res](https://github.com/HanRuize/DSAR-CF/blob/master/figs/res.png)

In horizontal-view videos, it is common to have subjects with mutual occlusion and being out-of-view. In this case, existing online trackers, e.g., DMAN can not associate the long-term lost subjects when they reappear in the view. Two examples are shown in Fig. 6.
The top two rows show the case of mutual occlusions. From the top view at frame #180, we can find that two subjects (ID
number 2, 3) are occluded by others and DMAN switches the ID of them when they reappear in the filed of view at frame #210. Our method keeps the original ID number. Similarly, we focus on the key subject (ID number 4) which goes out of view at frame #165 in the horizontal view. We can find that this subject is reassigned to a new ID number by DMAN. Our approach gets the original ID number of the target, which is consistent to its ID number in the top view


Dataset: https://pan.baidu.com/s/18TEX2OarO6KwikoWZM-GyA

Code: Coming soon...
