通常bioconductor和cran的R包都能装上。

github的有时候会有点难装，直接的`devtools::install_github("sqjin/CellChat")`直接报错。

通常都是依赖的包编译问题，这时就必须`devtools::install_github("sqjin/CellChat", dependencies = F, quiet = F)`

把依赖的包通过bioconductor和cran装上，最后就可以直接装成功了。

最难的是一个系统层面的问题，比如rjava的安装，得把编译层面的都装好才行。

