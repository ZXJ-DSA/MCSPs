## Introduction

This is the implementation code of the TKDE paper "I/O-Efficient Multi-Criteria Shortest Paths Query Processing on Large Graphs". Please refer to the paper for the algorithm details.

## Algorithms

1.` mcspem.cpp` is then entrance of all the proposed algorithms (i.e. OHP, MHP, BMHP, BMHPS).
2. Since STXXL vector requires const values for specify its parameters, you may need to manually specify the memory size each algorithm uses to store graph data in `head.h` if any. Sorry for the inconvenience.


## Data

#### Datasets

The original datasets can be downloaded from `http://users.diag.uniroma1.it/challenge9/download.shtml`, `http://konect.cc/` and `https://manliodedomenico.com/data.php`. Note that preprocessing is indispensable before running. 
For your convenience, the example processed graph data (US) are provided in following link for your reference. Link: `https://pan.baidu.com/s/1r2hVuA4K7r97rkYtsS57OQ` Passcode: `b5wm`

## Dependency

`g++` and `STXXL 1.4.1`.


## Reference

If you found any of our implimentations useful, please cite the appropriate references, listed below:

```
@article{zhou2024efficient,
  title={I/O-Efficient Multi-Criteria Shortest Paths Query Processing on Large Graphs},
  author={Zhou, Xinjie and Huang, Kai and Li, Lei and Zhang, Mengxuan and Zhou, Xiaofang},
  journal={IEEE Transactions on Knowledge and Data Engineering},
  year={2024},
  publisher={IEEE}
}
```
