# SC-LeGO-LOAM (noted)
## NEWS 
-（2026.3.17） 第一次学习scan context代码，为Scancontext.cpp加了注释。


## Features 
- Light-weight: a single header and cpp file named "Scancontext.h" and "Scancontext.cpp"
    - Our module has KDtree and we used <a href="https://github.com/jlblancoc/nanoflann"> nanoflann</a>. nanoflann is an also single-header-program and that file is in our directory.
- Easy to use: A user just remembers and uses only two API functions; `makeAndSaveScancontextAndKeys` and `detectLoopClosureID`.
- Fast: The loop detector runs at 10-15Hz (for 20 x 60 size, 10 candidates)


## Scan Context integration

- For implementation details, see the `mapOptmization.cpp`; all other files are same as the original LeGO-LOAM.
- Some detail comments
    - We use non-conservative threshold for Scan Context's nearest distance, so expect to maximise true-positive loop factors, while the number of false-positive increases.
    - To prevent the wrong map correction, we used Cauchy (but DCS can be used) kernel for loop factor. See `mapOptmization.cpp` for details. (the original LeGO-LOAM used non-robust kernel). We found that Cauchy is emprically enough.
    - We use both two-type of loop factor additions (i.e., radius search (RS)-based as already implemented in the original LeGO-LOAM and Scan context (SC)-based global revisit detection). See `mapOptmization.cpp` for details. SC is good for correcting large drifts and RS is good for fine-stitching.
    - Originally, Scan Context supports reverse-loop closure (i.e., revisit a place in a reversed direction) and examples in <a href="https://github.com/kissb2/PyICP-SLAM"> here (py-icp slam) </a>. Our Scancontext.cpp module contains this feature. However, we did not use this for closing a loop in this repository because we found PCL's ICP with non-eye initial is brittle. 

## How to use 
- Place the directory `SC-LeGO-LOAM` under user catkin work space 
- For example, 
    ```
    cd ~/catkin_ws/src
    git clone https://github.com/OlivaYao/SC-LeGO-LOAM-noted-.git
    cd ..
    catkin_make
    source devel/setup.bash
    roslaunch lego_loam run.launch
    ```

## MulRan dataset 
- If you want to reproduce the results as the above video, you can download the <a href="https://sites.google.com/view/mulran-pr/home"> MulRan dataset </a> and use the <a href="https://sites.google.com/view/mulran-pr/tool"> ROS topic publishing tool </a>.   


## Dependencies
- All dependencies are same as LeGO-LOAM (i.e., ROS, PCL, and GTSAM).
- We used C++14 to use std::make_unique in Scancontext.cpp but you can use C++11 with slightly modifying only that part.
