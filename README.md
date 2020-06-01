#  CFD with Fortran
**此项目存储所选课程《非线性方程组数值解法》中我所编写的作业题**  
*因为是作业，所以程序都是只有主程序，没有做继续的细化，也没有把用到的模块单独列出来*
## 1D Burgers' equ
首先是对于一维Burgers' equ的简单求解，对于连续解情况，没有加入限制器，最多达到O(Δx^3)  
![Buergers'equ]()
对于最基本的限制器，这里采用了TVD与TVB进行编写  

## WENO limiter
对于限制器这里我只提供了WENO，其他的限制器用的比较少，这里就不赘述了
对于处理方法分别采用了**有限体积法（FV method）**以及**有限差分法（FD method）**  

## 1D Euler 方程组的求解
*这是最后的大作业，分两问，第一问是连续解情况，第二问是间断解情况*  
<img scr="https://github.com/lmd01/CFD_with_Fortran/blob/master/Euler_equations/6C3DED3299C5FD79DBE2DDA14BA9E3A6.JPG" width="200" height="200" alt="Euler"/><br/>  
这里给出用Mathematica软件画出的间断解结果图，间断解为初值发展0.6秒后的结果  
![image1](https://github.com/lmd01/CFD_with_Fortran/blob/master/Euler_equations/%E9%97%B4%E6%96%AD%E8%A7%A3%E7%A4%BA%E6%84%8F%E5%9B%BE.jpg)
