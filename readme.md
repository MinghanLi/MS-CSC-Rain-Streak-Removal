README

[Video Rain Streak Removal By Multiscale Convolutional Sparse Coding](https://openaccess.thecvf.com/content_cvpr_2018/papers/Li_Video_Rain_Streak_CVPR_2018_paper.pdf)
[project](https://sites.google.com/view/cvpr-anonymity) [dataset](http://gr.xjtu.edu.cn/web/dymeng/2)

================================================================================================
Version 1.0, 28-Apr-2018

The code has been tested with MATLAB 2017b on Windows.

Use of this code is free for research purposes only. Shall not be used for commercial purposes! 
All copyrights belong to the original authors. The technology has applied for patents. If you 
want to purchase the patents for commercial purposes, please contact the corresponding author: 
Deyu Meng, dymeng@mail.xjtu.edu.cn. Thank you!


# Reference: 

Minghan Li, Qi Xie, Qian Zhao, Wei Wei, Shuhang Gu, Jing Tao, Deyu Meng*, Video Rain Streak Removal By Multiscale Convolutional Sparse Coding. CVPR, 2018.


# Usage: 

1. Unpack the contents of the compressed file to a new directory.

2. Run the DemoY.


# Dataset: 

Synthetic rain video: 

'park.mat', 'highway.mat', 'man.mat', corresponding to 
Fig. 3-5 in the reference. 

Real rain video: 'wall.mat', 'compfinal.mat', 'night.mat' corresponding to Fig. 6-8 in the reference, respectively.

 
More experiments on real videos are shown in [project](https://sites.google.com/view/cvpr-anonymity) .
The dataset can be downloaded from [here](http://gr.xjtu.edu.cn/web/dymeng/2) .



# Acknowledgments 

Code partly borrows from:

[1] Wei Wei, Liyuan Yi, Qi Xie, Qian Zhao, Deyu Meng, Zongben Xu, Should We Encode Rain Streaks 
in Video as Deterministic or Stochastic? ICCV, 2017.

[2] Shuhang Gu, Deyu Meng, Wangmeng Zuo, Lei Zhang. Joint Convolutional Analysis and Synthesis Sparse Representation for Single Image Layer Separation. ICCV, 2017.
Thanks for sharing!

# Our latest paper
We propose a new online rain/snow removal method from surveillance videos by fully encoding the dynamic statistics of both rain/snow and background scenes in a video along time into the model, and realizing it with an online mode to make it potentially available to handle constantly coming streaming video sequence. 

[Minghan Li](https://scholar.google.com/citations?user=LhdBgMAAAAAJ&hl=en&oi=ao),
[Xiangyong Cao](https://scholar.google.com/citations?user=IePM9RsAAAAJ&hl=en),
[Qian Zhao](https://scholar.google.com/citations?user=vM6yGTEAAAAJ&hl=en),
[Lei Zhang](https://scholar.google.com/citations?user=tAK5l1IAAAAJ&hl=en&oi=ao),
[Deyu Meng](https://scholar.google.com/citations?user=an6w-64AAAAJ&hl=en&oi=ao)

Please go to the [Homepage](https://github.com/MinghanLi/OTMSCSC_matlab_2020) to obtain more information about our work.

```
@inproceedings{li2018video,
  title={Video rain streak removal by multiscale convolutional sparse coding},
  author={Li, Minghan and Xie, Qi and Zhao, Qian and Wei, Wei and Gu, Shuhang and Tao, Jing and Meng, Deyu},
  booktitle={Proceedings of the IEEE conference on computer vision and pattern recognition},
  pages={6644--6653},
  year={2018}
}
@article{Li2021OnlineRR,
  title={Online Rain/Snow Removal From Surveillance Videos},
  author={Minghan Li and Xiangyong Cao and Q. Zhao and L. Zhang and Deyu Meng},
  journal={IEEE Transactions on Image Processing},
  year={2021},
  volume={30},
  pages={2029-2044}
}
```

# Contact

If you have any questions, please feel free to contact Minghan Li(liminghan0330@gmail.com).
