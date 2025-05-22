# 做peanut的pySCENIC分析的故事

写在前面：空洞的去记录一个工作和流程很容易遗忘，将具体案例/项目记下来，后面再看就会深刻许多，也好懂很多。该工作的难点如下：pySCENIC环境获得；可视化用的R版本的SCENIC环境获得；最大的难点应该就是peanut作为非模式物种缺乏输入的背景文件（tf列表`.txt`，motif信息`.feather`，tf_motif对应信息`.tbl`）。

环境搭建
pySCENIC环境：
SCENIC环境：
blastp序列比对环境：
motif处理环境：

过程：
下载拟南芥的TF列表和对应的蛋白质序列文件；
找到peanut蛋白质序列文件；
使用blastp将peanut蛋白质序列比对到拟南芥上，拿到对应关系；
找到peanut单细胞数据建index基因组的fasta和gtf；
下载拟南芥motif的meme文件，tf_motif对应信息的txt文件；
对tf_meme对应信息的txt文件进行简单处理拿到拟南芥的tbl文件；
拟南芥meme文件提取motif信息，拿到各个motif特征文件；
对peanut的genome文件处理，提取启动子序列保存未fasta文件；
利用[create_cistarget_motif_databases.py](https://github.com/aertslab/create_cisTarget_databases/blob/master/create_cistarget_motif_databases.py)将提取到的peanut启动子序列比对到motif信息上，拿到peanut存在的motif信息，即ranking.feather文件；
基于peanut和拟南芥比对信息，将拟南芥的tf_list和tbl信息进行替换，没匹配的拟南芥信息删除。
