# Adjustment-Of-Traverse-Network
Matlab导线网自动间接平差

- 导线网自动间接平差
- 配套于武汉大学测绘学院测量平差仿真辅助教学系统V2.0
- 仅限于平差课程学习与交流
- Email: lijintao.alex@qq.com

## 三种程序使用方式

### 下载release包安装为Matlab App

[Measurement Ajustment for Traverse Network.mlappinstall](https://github.com/AkexStar/Adjustment-Of-Traverse-Network/releases)

### 下载m文件直接运行

[MeasurementAjustmentforTraverseNetwork.m](https://github.com/AkexStar/Adjustment-Of-Traverse-Network/blob/main/MeasurementAjustmentforTraverseNetwork.m)

### clone文件夹到本地用Matlab App Designer打开运行

[main.mlapp](https://github.com/AkexStar/Adjustment-Of-Traverse-Network/tree/main/AppDesignerCode)

## 运行效果截图

![screenshot](https://user-images.githubusercontent.com/55226358/205250148-a5a1e11a-8a2d-452a-a038-d036d27fde84.png)

![image](https://user-images.githubusercontent.com/55226358/205250351-046d6c99-f0e4-4573-ac36-a8f6a26f6029.png)

## 程序使用说明

- 本程序配套于武汉大学测绘学院测量平差仿真辅助教学系统V2.0
- 本程序中“边数据”指观测边边长数据
- 本程序中“角数据”指观测角角度数据

### 使用流程如下：

1. 将角数据与边数据粘贴进指定文本框
2. 点击读入角数据与边数据按钮
3. 若数据读入成功 提示灯变绿且自动显示2个已知点点号
4. 填入已知点坐标数据
5. 点击平差计算按钮

**边数据请粘贴进上方指定区域内**

内容包含边对应点号与边长(m),下面为样例数据格式示例：

0 -- 1	97.711

1 -- 2	69.81

2 -- 3	78.536

3 -- 4	72.12

4 -- 5	67.714

5 -- 6	已知值

6 -- 0	65.505

3 -- 7	77.436

7 -- 8	40.1

8 -- 6	68.199

**角数据请粘贴进上方指定区域内**

内容包含对应角号与角度(°′″),下面为样例数据格式示例:

∠0 1 2	236°00′33.5″

∠1 0 6	130°33′18.9″

∠1 2 3	225°51′8.5″

∠2 3 4	231°45′21″

∠2 3 7	281°25′11.9″

∠3 4 5	213°46′6.1″

∠4 3 7	49°40′8.4″

∠4 5 6	246°24′48″

∠5 6 0	236°45′44.4″

∠5 6 8	323°59′33.9″

∠0 6 8	87°13′30.7″

∠3 7 8	174°53′58.8″

∠7 8 6	199°35′58.9″

**已知点数据需手动填入，下面为样例数据：**

5 x=164.668(m), y=112.313(m)

6 x=274.722(m), y=136.706(m)

本程序默认测边误差：1/2000

默认测角误差：12″
