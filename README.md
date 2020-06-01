#  CFD with Fortran
**此项目存储所选课程《非线性方程组数值解法》中我所编写的作业题**  
*因为是作业，所以程序都是只有主程序，没有做继续的细化，也没有把用到的模块单独列出来*
## 1D Burgers' equ
首先是对于一维Burgers' equ的简单求解，对于连续解情况，没有加入限制器，最多达到O(Δx^3)  
对于最基本的限制器，这里采用了TVD与TVB进行编写  

## WENO limiter
对于限制器这里我只提供了WENO，其他的限制器用的比较少，这里就不赘述了
对于处理方法分别采用了**有限体积法（FV method）**以及**有限差分法（FD method）**  

## 1D Euler 方程组的求解
*这是最后的大作业，分两问，第一问是连续解情况，第二问是间断解情况*
这里给出用Mathematica软件画出的间断解结果图，间断解为初值发展0.6秒后的结果  
![image](https://github.com/lmd01/CFD_with_Fortran/Euler equations/间断解示意图.jpg)
