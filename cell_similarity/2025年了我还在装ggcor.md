# 2025年了我还在装ggcor

```shell
wget https://github.com/Github-Yilei/ggcor/archive/refs/heads/main.zip
unzip master-main.zip
Rscript -e 'devtools::install_local("ggcor-main")'
Rscript -e 'library(ggcor)'
Rscript -e "devtools::install_github('thomasp85/ambient')"
Rscript -e 'library(ambient)'
```
Thanks to them:

[共享服务器上ggcor包怎么装](https://mp.weixin.qq.com/s/-3RSe5TMwjRk-2BGCBRTtg)
[利用ggcor包绘制相关性组合图及环状热图](https://www.jianshu.com/p/c5f19a166bc8)
[ggcor包作图 | 相关性热图 | mental分析图](https://zhuanlan.zhihu.com/p/507384776)
[ggcor 的环形热图](https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ)

这个事情我干了一天……