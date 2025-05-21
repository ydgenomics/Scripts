[同源基因比对：顶刊都在用的跨物种分析](https://mp.weixin.qq.com/s/jv2Z8NVWZwzeVjmm5c9NNg)
[使用ensemble得到两个物种的同源基因：从BioMart中下载拟南芥和小麦的直系同源基因（Orthologues）](https://www.jianshu.com/p/5de2c98797f2)

同源基因和直系同源基因的区别？

学习使用blast做两个不同物种的蛋白质序列找同源基因
BLASTP本身主要用于单序列与数据库的比对，而不是用于多序列比对
[史上最详细的blast安装附视频](https://mp.weixin.qq.com/s/rEBqjN-fGOp_loTmyEuMJA)
```shell
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
# environment location: /software/ncbi-blast-2.16.0+/bin
vim ~/.bashrc
export PATH=/software/ncbi-blast-2.16.0+/bin:$PATH
source ~/.bashrc
blastp -h
```

blastp不适合多序列比对，使用Clustal加快序列比对
[Clustal Omega—广泛使用的多序列比对工具](https://mp.weixin.qq.com/s/f9pEFWJJoNCqlFEfd77aOA)