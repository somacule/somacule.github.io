---
layout: post
title: Get all atoms RDF
date: 2020-03-07 20:20:23 +0900
category: calculation
---
# RDF定义
> [径向分布函数](https://baike.baidu.com/item/%E5%BE%84%E5%90%91%E5%88%86%E5%B8%83%E5%87%BD%E6%95%B0/12723225?fr=aladdin)（Radial distribution function）通常指的是给定某个粒子的坐标，其他粒子在空间的分布几率（离给定粒子多远）。所以径向分布函数既可以用来研究物质的有序性，也可以用来描述电子的相关性。
>
> 如果给定粒子当做原点，体系平均粒子数密度为*r<sub>ho</sub>*=N/V，则距原点为r处的局部时间平均的密度为*r<sub>ho</sub>*×g(r) 。这是对均匀的各向同性系统的简化定义。
>
> 简言之，这是对于距参考粒子距离为*r*处找到粒子的相对概率的测量，参考态是理想气体。一般的算法是计算在距参考原子*r*到*r*+d*r*这样的壳层里有多少粒子。如下图，深红为参考粒子，蓝色为找到的粒子，在*r*到*r*+d*r*的范围（虚线表示）。
>
> 通常径向分布函数RDF表示为g(*r*)

<span><div style="text-align: center;">

![中心原子的r到dr之间](https://ss0.bdstatic.com/70cFvHSh_Q1YnxGkpoWK1HF6hhy/it/u=2380327841,484091491&fm=26&gp=0.jpg)

</div></span>


# RDF计算

RDF可以从多种理论和实验方法得到，包括但不限于经典/从头计算分子动力学（MD/AIMD），同步辐射X射线拓展边吸收谱（EXAFS）等。其中，经典分子动力学由于计算量较小，精度尚可，且模拟体系较大，可达成千上万个原子，因此常用模拟平衡后的体系，进行RDF计算。常见的计算软件自带的计算方法通常只能得出一类原子的RDF数据，统计平均化了各个原子的局域配位环境。本脚本可以通过使用gromacs内建RDF程序导出单个原子的RDF，最后直接求解单个原子的配位环境。

# Code

```python
#! /usr/bin/python
import os
import sys
import numpy as np
import scipy
from scipy import integrate
global listofmd
listofmd = ['4','100','100','0.330','0.300']
#Box-r,Cl-num,O-num,Cl-cutoff-r,O-cutoff-r
global coord_Cl, coord_O
coord_Cl = []
coord_O = []
#-----------------------------------------------------------
def intro(entry):
	global listofmd
	global coord_Cl, coord_O
	f_org = open(entry)
	entry_n = entry.replace('.xvg','')
	lines = f_org.readlines()
	con1 = False
	conCl = False
	conO = False
	str1 = []
	str2 = []
	str3 = []
	int_str2 = []
	int_str3 = []
	for line in lines:
		if line.find('0.000') != -1:
			con1 = True
		if con1:
			list1 = line.split()
			str1.append(float(list1[0]))
			str2.append(cal_Cl(list1[0],list1[1]))
			str3.append(cal_O(list1[0],list1[2]))
			if list1[0]==listofmd[3]:
				conCl = True
			if list1[0]==listofmd[4]:
				conO = True
		if conCl:
			coord_Cl.append(scipy.integrate.simps(str2, str1, even='avg'))
			conCl = False
		if conO:
			coord_O.append(scipy.integrate.simps(str3, str1, even='avg'))
			conO = False
		if conCl and conO:
			break
	f_org.close()

#calculating r^2*g(r)*4pi/R^3 from str to float
def cal_Cl(r,gr):
	global listofmd
	r2gr = (float(listofmd[1])*4*3.14159*float(r)**2*float(gr))/(float(listofmd[0])**3)
	return r2gr
def cal_O(r,gr):
	global listofmd
	r2gr = (float(listofmd[2])*4*3.14159*float(r)**2*float(gr))/(float(listofmd[0])**3)
	return r2gr

def write_ClO():
	global coord_Cl, coord_O
	f_new = open('coordinated_Cl_O.csv','w')
	for i in range(0,len(coord_Cl)):
		f_new.write(str(coord_Cl[i])+','+str(coord_O[i])+'\n')
	f_new.close()              



list = os.listdir(os.curdir)
i = 0
while True:
	if i >= len(list): break
	if not list[i].endswith('.xvg'):
		list.pop(i)
		continue
	i +=1
	
for f in list:
	intro(f)
write_ClO()
```
