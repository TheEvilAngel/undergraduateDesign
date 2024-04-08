最大的bug
1. 共轭没弄好->ifft2后是复数，结果为了方便直接real了（可以改进但懒了）
2. 论文写错了，age实际上是-0.75次方【所以有风域画图注释的纠结】，其他的数值没变大差不差也没改^ ^
3. 论文写的ifft2前的系数应该和海面面积联合起来->功率守恒，我当时弄完发现两个画图差不多，也没改^ ^
PS：2~3可以试着改改，1不用弄了

generateSeaSurface——不含时间因子，最古老版本，绝对没问题【推荐这个】
generateSeaSurface_time_new——k是w的函数
generateSeaSurface_ship——画视频用的
eg——知乎学习海面代码用的