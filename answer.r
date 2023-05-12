
###################################截取中间值
dat1 <- c('human(display_long)|uniprotkb:ESR1(gene name)')
dat2 <- c('human(display_long)|uniprotkb:TP53(gene name)')
dat3 <- c('human(display_long)|uniprotkb:GPX4(gene name)')
dat4 <- c('human(display_long)|uniprotkb:ALOX15(gene name)')
dat5 <- c('human(display_long)|uniprotkb:PGR(gene name)')

dat <- rbind(dat1,dat2,dat3,dat4,dat5)
a <- substring(dat,31,) #截去前面
a <- substr(a,1,nchar(a)-11)  #截去后面


####################################定义域问题
rm(list=ls())
library(ISLR2)
u=-1.5
D=1 #u小于等于正负1
D=0 #u大于正负1
M=D*(15/16)*(1-u^2)^2

qxart1=function(x){
  if(abs(x)<=1){D=1}else{D=0}
  M=D*(15/16)*(1-x^2)^2
  return(M)
}



x=0
1-x^2

plot(quart1,
     xlim=c(-1,1),
     type="p")

rule=fuction(u){u>-1&&u<1}
a=c(-1,1)


quart1=function(x){
  ifelse(abs(x)<=1,M=(15/16)*(1-x^2)^2,M=0)
  return(M)
  }







###########################geo下载不了数据
options(stringsAsFactors = F)
f='GSE31684_eSet.Rdata'
library(BiocGenerics)
library(Biobase)
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE31684', destdir=".",
                 AnnotGPL = F, ## 注释文件
                 getGPL = F,GSEMatrix = T,parseCharacteristics = T) ## 平台文件
  save(gset,file=f) ## 保存到本地
}
load('GSE31684_eSet.Rdata')



########################图例无法显示问题
#输入数据
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
opar <- par(no.readonly=TRUE)

#绘制图形
plot(dose, drugA, type="b",
     pch=15, lty=1, col="red", ylim=c(0, 60),xlim = NULL,
     main="Drug A vs. Drug B",
     xlab="Drug Dosage", ylab="Drug Response",)
lines(dose, drugB, type="b",
      pch=17, lty=2, col="blue")
abline(h=c(30), lwd=1.5, lty=2, col="gray")

#添加图例
legend("topleft", inset=.05, title="Drug Type", c("A","B"),lty=c(1, 2), pch=c(15, 17), col=c("red", "blue"))


##############
a <- read.table("GSE188952.txt",sep = "\t",header = F,row.names = 1)
a <- t(a)
b <- a[,1]
write(b,file = "gene.txt",sep = "\t")

geneFile="gene_happy.txt"        #基因列表文件
gene=read.table(geneFile,sep="\t",header=F,check.names=F)
out=gene[,1]
write(out,file = "out.txt",sep = "\t")

###############
aa=1111
Fre=2222
file4=as.character(paste(aa,"/",aa,"-",Fre,"-gene",".txt", sep=""))

###############
library(car)
library(ggplot2)
state <- as.data.frame(state.x77[,c('Murder','Population','Frost','Income')])
fit <- lm(Murder~Population+Frost+Income,data = state)
qqPlot(fit,
       row.names(state),
       simulate=T,col="black", col.lines="red",main="qq plot",add.line=FALSE)


###############
x <- 1:10
y <- x^2
ci_l <- x^2 - 0.5 * x
ci_r <- x^2 + 0.5 * x

dat_plot <- data.frame(x, y, ci_l, ci_r)
ggplot(dat_plot, aes(x = x)) +                  # x轴在此处添加，目的为了置信区间与拟合线共享同一个x
  geom_ribbon(aes(ymin = ci_l, ymax = ci_r)) +  # 添加置信区间
  geom_line(aes(y = y))                         # 添加拟合线

set.seed(357) 
dat <- rgamma(20, shape=2, scale=2) 
dev.new()
qqPlot(dat)

################
equation <- function(a, b, k){
  b <- rev(b)
  concat <- numeric(length(a))
  for (i in seq_along(a)){
    concat[i] <- paste0(a[i],b[i])
  }
  test <- concat < k
  result <- sum(test == TRUE)
  result
}
##a取三位数
a<- c(1, 2, 3)
b<- c(1, 2, 3)
k<- 31
equation(a, b, k)


###################y轴0和100倒置，换成百分数
library(ggplot2)
head(iris)
#绘制变量Sepal.Length的累积分布
ggplot(iris,aes(x = Sepal.Length))+
  stat_ecdf(color = "red")+
  scale_y_reverse(labels = scales::percent) #y轴上下倒置+转换百分比


###################设置y轴均匀刻度
library(ggplot2)
head(iris)
#绘制变量Sepal.Length的累积分布
ggplot(iris,aes(x = Sepal.Length))+
  stat_ecdf(color = "red")+
  ylim(0,15)


###################
a <- iris
length(iris)

##################
library(survival)
library(rms)
data <- colon
dd=datadist(data)
options(datadist="dd")
f2 <- psm(Surv(time, status) ~ surg+node4+obstruct, data, x=T, y=T, dist='lognormal',time.inc = 180)
f2
cal1 <- calibrate(f2,
                  cmethod='KM',
                  method="boot",
                  u=180, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=61, #每次抽样的样本量，
                  B=100) #抽样次数

plot(cal1,lwd=2,lty=1,
     conf.int=T,# 是否显示置信区间
     errbar.col="black",#直线曲线bar颜色
     col="black", # 曲线颜色
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted Probability",
     ylab="Actual Probability",
     subtitles = F)#不显示副标题
#abline(0,1,col="black",lty=1,lwd=1)
lines(cal1[,c("predy","calibrated.orig")],type="l",lwd=1,col="red",pch=16)
lines(cal1[,c("predy","calibrated.corrected")],type="l",lwd=1,col="blue",pch=16)



#################
install.packages("xlsx")
library(xlsx)
a <- iris

write.xlsx(shiyan_1_1,"实验8+学号+姓名+R语言数据匹配检查表-学生表.xlsx")

#################提取疾病共同基因
library(dplyr)
dat <- read.table("test.txt",sep = "\t",header = T)
dat
a <- dat[dat$disease=="a",]
b <- dat[dat$disease=="b",]
c <- dat[dat$disease=="c",]

ab <- inner_join(a,b,by="gene")
bc <- inner_join(b,c,by="gene")
ac <- inner_join(a,c,by="gene")


################取n-1项的均值
a <- (1:10)
a
mean(a[1:(length(a)-1)])

###############
install.packages("remotes")
remotes::install_github("audrey-b/BUGSnet")

##############
dat <- read.table("deg.txt",header = 1)
by(dat$a,dat$klg1,shapiro.test) 

##############
for (i in 1:4898) {
  if(wine[i,12]>6)cha[i]="good"
  else if(wine[i,12]>5)cha[i]="mid"
  else cha[i]="bad"}

iris <- iris
iris1 <- iris[3,3]

#############glm逻辑回归两个报错处理
library("ggplot2")
data<-iris[1:100,]
samp<-sample(100,80)
names(data)<-c('sl','sw','pl','pw','species')
testdata<-data[samp,]
traindata<-data[-samp,]
lgst<-glm(testdata$species~pl,binomial(link='logit'),data=testdata)
summary(lgst)

#先计算测试集
lgst<-glm(testdata$species~pl,binomial(link='logit'),data=testdata,control=list(maxit=27))
testdata$p<-predict(lgst,type='response')

qplot(seq(-2,2,length=80),sort(p),col='predict')
#检查测试集pl值
testdata$y <- c(1:80)
qplot(pl,y,data =testdata,colour =factor(species))
testdata

############合并三个矩阵
BP <- matrix(1:12,ncol = 4,byrow = T)
CC <- matrix(1:8,ncol = 4,byrow = T)
MF <- matrix(1:4,ncol = 4,byrow = T)
GO<-rbind(BP,CC,MF)

############STRINGdb报错
#BiocManager::install("STRINGdb")
library(STRINGdb)
string_db <- STRINGdb$new( version="11.5", #数据库版本。截止2022.5.24最新为11.5
                           species=9606,   #人9606，小鼠10090 
                           score_threshold=400, #蛋白互作的得分 默认400, 低150，高700，极高900
                           input_directory="") #可自己导入数据

############
install.packages("remotes")
remotes::install_github("PelzKo/immunedeconv2")

############
data <- iris
p<-ggPLOT(data,AES(x=dAtA5sinstallmEntfill=grade))

###########把循环出来的数值变成变量
x<-1:100
x
for(i in x*10){
  y<-mean(rbinom(i,1,0.3))
  y
  print(y)
}

##########
version
pkgbuild::check_rtools(TRUE)
pkgbuild::has_rtools(TRUE)
pkgbuild::has_build_tools(debug=TRUE)
devtools::find_rtools(debug=TRUE)
library(devtools)
library(ggbiplot)
install_github("vqv/ggbiplot")

#########
devtools::install_github("NightingaleHealth/ggforestplot")
rm (list = ls ())
library(meta)
library(forestplot)
library(ggplot2)
library(tidyverse)
library(ggforestplot)
data(Olkin95)
m1 <- metabin(event.e, n.e, event.c, n.c,data = Olkin95, subset = c(41, 47, 51, 59),sm = "RR", method = "I",studlab = paste(author, year))
forest(m1);m1
dev.new()

#########数据操作和绘图
data <- read.table("cgss2010",sep = "\t",header = T,row.names = 1)
tar <- c("ID",
         "S5",
         "S41",
         "A2",
         "A3A",
         "A3B",
         "A3C",
         "A4",
         "A7A",
         "A15",
         "A8A",
         "A36",
         "L6A")
sample(data, 50, replace = FALSE, prob = NULL)
set.seed(2019101027)
data1 <- data[sample(nrow(data), 50), ]

########取交集
drugA
drugB
intersect <- Reduce(intersect,list(drugA,drugB))
intersect

########添加分组
library(stringr)
data <- iris
data
a <- iris[1:5,]
b <- iris[51:55,]
c <- iris[101:105,]
dat1 <- rbind(a,b,c)
x <- ifelse(str_detect(dat$Species,"setosa"),"group1","group2")
dat1$group <- x
dat2 <- dat1
dat2
library(tidyverse)
length(which(dat2$Species=="setosa")) #每个分组的个数
a1 <- mutate(dat1[str_detect(dat2$Species,"setosa"),],group=rep("group1",5))
b1 <- mutate(dat1[str_detect(dat2$Species,"versicolor"),],group=rep("group2",5))
c1 <- mutate(dat1[str_detect(dat2$Species,"virginica"),],group=rep("group3",5))
x <- rbind(a1,b1,c1)

#######随机森林
library("party")
library("randomForest")
# Create the forest.
output.forest <- randomForest(nativeSpeaker ~ ., 
                              data = readingSkills)

# View the forest results.
print(output.forest) 

# Importance of each predictor.
print(importance(output.forest,type = 2))

#########table1
library(table1)
library(boot)
dat <- read.csv("bird.csv",row.names = 1)
colnames(dat)
table1(~ Num_1990 + ASDR_1990 + Num_2019 + ASDR_2019 + EAPC_CI  | location, data=dat)

dat1 <- read.csv("bird2.0.csv")
colnames(dat1)
a <- table1(~ China+Global+High.middle.SDI+High.SDI+Low.middle.SDI+Low.SDI+Middle.SDI  | year*location, data=dat1)

#########
plot(1:30)
plot(1:30,main="图形元素设置演示")

library(ggplot2)
a <- plot(1:30)
a+labs(title = "Title of the plot",
     subtitle = "Subtitle of the plot",
     caption = "This is the caption",
     tag = "Fig. 1")+  
  theme(plot.title = element_text(family = "serif", #字体
                                   face = "bold",     #字体加粗
                                   color = "red",      #字体颜色
                                   size = 15,          #字体大小
                                   hjust = 0.5,          #字体左右的位置
                                   vjust = 0.5,          #字体上下的高度
                                   angle = 0,          #字体倾斜的角度
))

text(30,1,"year",pos=4)
plot(1:30, ann = F, xaxt = "n", yaxt ="n")
axis(1, -4:4, -4:4, mgp = c(3, 2, 0)) 
axis(1, -4:4, 1, padj = -1)

##########
a <- iris[1:75,]
a <- a[,4:5]
b <- iris[76:150,]
b <- b[,3:5]

colnames(a)
c <- merge(a,b,by='Species',all = F)
c <- merge(a,b,by=intersect(names(a)[2],names(b)[3]))
c <- merge(a,b ,by='Species',all.y=TRUE,sort=TRUE) 

########gg.gap截断y轴
# install.packages("gg.gap")
library(gg.gap)
library(ggplot2)
data <-
  data.frame(x = c("Alpha", "Bravo", "Charlie", "Delta"),
             y = c(200, 20, 10, 15))
#画图
p1 = ggplot(data, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = 'identity', position = position_dodge(),show.legend = FALSE) +
  theme_bw() +
  labs(x = NULL, y = NULL)

p1 = p1 +  theme_classic()
p1
p2 = gg.gap(plot = p1,
            segments = c(25, 190),
            tick_width = 10,
            rel_heights = c(0.25, 0, 0.1),# 设置分隔为的三个部分的宽度
            ylim = c(0, 200)
)

p2

#######
library(sjPlot)
library(tinytex)
library(compareGroups)

af=read.csv("bird.csv", header=T, check.names=F,fill = T)
af$`KEGG enrichment`=paste0(af$id,": ",af$KEGG)
af=af[,-2]
tab_df(a,file="sxb9.html")

#######
data = diamonds
g <- ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, y = stat(prop), group = x, fill = cut))
g
ggplot(diamonds)+geom_bar(aes(x=clarity, fill=cut))
p4 <- g + geom_bar(aes(fill=class)) +
  scale_fill_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", 
                               "#c7eae5", "#5ab4ac"))
p4

########
# Load the party package. It will automatically load other
# required packages.
library(party)
# Print some records from data set readingSkills.
print(head(readingSkills))
library(party)
library(randomForest)
# Create the forest.
output.forest <- randomForest(nativeSpeaker~.,data = readingSkills)


########
library(pheatmap)
data <- readingSkills[,-1]
data=log2(data+1)
pdf("heatmap.pdf", width=7.5, height=5)
pheatmap(data,
#annotation=Type,
annotation_colors = colorList,
color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
cluster_cols =F,
cluster_rows =T,
scale="row",
show_colnames=F,
show_rownames=T,
fontsize=6,
fontsize_row=6,
fontsize_col=6)


########使用R调用腾讯地图API批量获取地点经纬度坐标
install.packages('jsonlite')
library(jsonlite)
key="KV3BZ-JNICX-ZN24F-75PQ7-CGUET-EXFU6"

#地址补全
ad="佛山市白坭中学"
c <- fromJSON(paste0("https://apis.map.qq.com/ws/smart_address/address_complete?address=",ad,"&key=",key))
resc <- c(c$result$completed_address,
          paste0(c$result$province$name,
                 c$result$city$name,
                 c$result$district$name,
                 c$result$town$name,
                 c$result$road$name,
                 c$result$village$name))
resc #标准地址，具体地址

#分析地标建筑
ad2="佛山市南海中学"
d <- fromJSON(paste0("https://apis.map.qq.com/ws/smart_address/place_analy?address=",ad2,"&key=",key))
resd <- c(d$result$pois$address,d$result$pois$category)
print(resd)
sub1 <- subset(as.data.frame(d$result$sub_pois), select=c(address,title,category,location))  #矩
#sub1 <- paste0(d$result$sub_pois$address,d$result$sub_pois$title,d$result$sub_pois$category,d$result$sub_pois$location) #行
print(sub1)
view(sub1)

#ip定位，只能准确到区
ip="119.131.17.224"
b <- fromJSON(paste0("https://apis.map.qq.com/ws/location/v1/ip?ip=",ip,"&key=",key))
resb <- c(b$result$location$lng,b$result$location$lat)
print(resb)

#经纬度
la="广东省广州市番禺区大学城外环东路382号"
a <- fromJSON(paste0("https://apis.map.qq.com/ws/geocoder/v1/?address=",la,"&key=",key))
resa <- c(a$result$location$lng,a$result$location$lat)
print(resa)

#批量经纬度
la2<-c("广东省广州市番禺区大学城外环东路382号",
       "广州市番禺区广州大学城外环东路232号",
       "广东省广州市番禺区大学城外环东路132号",
       "广东省广州市番禺区大学城外环东路280号")
k <- list()
for (i in la2) {
  a <- fromJSON(paste0("https://apis.map.qq.com/ws/geocoder/v1/?address=",i,"&key=",key))
  resa2 <- c(a$result$location$lng,a$result$location$lat)
  k[[i]] <- print(resa2)
}
resa3 <- t(as.data.frame(k))
print(resa3)

#可视化
uni_geo <- get_geo_position(la2)
uni_result <- remapB(markPointData = data.frame(resa3),
                     markPointTheme = markPointControl(symbol = 'circle',
                                                       effect = T,
                                                       symbolSize = 12,
                                                       color = 'red'),
                     geoData = uni_geo
)

#############
input <- mtcars[,c('mpg','cyl')]
print(head(input))

boxplot(mpg ~ cyl, data = mtcars, xlab = "good",
        ylab = "bad", main = "abc")

#############
install.packages("ggtern")
library(ggtern)
data(USDA)
draw_key_crosshair_tern(data, params, size)
data(WhiteCells)
ggtern(WhiteCells,aes(G,L,M)) + 
  geom_density_tern(aes(color=Experiment)) +
  geom_point(aes(shape=Experiment)) +
  facet_wrap(~Experiment,nrow=2)

############
install.packages("plot3D")
install.packages("rgl")
library(plot3D)
library(rgl)
M = mesh(seq(-pi,pi,length = 50),seq(-pi,pi,length = 50))
u = M$x
v = M$y
r = 1
x = r*cos(u)*cos(v)
y = r*cos(u)*sin(v)
z = r*sin(u)
plot3d(x,y,z,col = 'blue',alpha = 0.3,type = 'l')
points3d(cos(pi/4)*cos(pi/4),cos(pi/4)*sin(pi/4),sin(pi/4),col = 'red',size = 5)

###########
library(rgl)
data(iris)
head(iris)

x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width
m <- pet.w <- iris$Petal.Width

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

open3d(windowRect = c(50, 50, 562, 562))# 打开一个新的 RGL device
rgl.points(x, y, z, color ="lightgray") # 散点图
rgl.quads(x, y, z, color ="lightgray")
rgl.bbox(color = "#333377") # 添加边界框装饰

# x, y, z : 对应点坐标的数值向量
#axis.col : 轴颜色
# xlab, ylab, zlab: 轴标签
# show.plane : 添加轴平面
# show.bbox : 添加边界框装饰
# bbox.col：边界框颜色。 第一种颜色是背景色； 第二种颜色是刻度线的颜色
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))

  rgl_add_axes <- function(x, y, z, axis.col = "grey",
                           xlab = "", ylab="", zlab="", show.plane = TRUE, 
                           show.bbox = FALSE, bbox.col = c("#333377","black"))
  { 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # 添加轴
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z); 
  rgl.lines(xlim, c(0, 0), c(0, 0),c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0),c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim,c(0, 0), color = axis.col)

  
  # 在每个轴的末端添加一个点以指定方向
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # 添加轴标签
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # 添加平面
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # 添加边界框装饰
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
  }

get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
cols <- get_colors(iris$Species, c("#999999", "#E69F00", "#56B4E9"))
rgl_init()
rgl.points(x, y, z, m, r = 0.2, color = cols) 
rgl_add_axes(x, y, z, m)
aspect3d(1,1,1)

rgl_init()
shapelist3d(tetrahedron3d(), x, y, z, size =  0.15, 
            color = get_colors(iris$Species))
rgl_add_axes(x, y, z, show.bbox = TRUE)
aspect3d(1,1,1)


library(plot3D)
library(rgl)
M = mesh(seq(-pi,pi,length = 50),seq(-pi,pi,length = 50))
u = M$x
v = M$y
r = 1
x = r*cos(u)*cos(v)
y = r*cos(u)*sin(v)
z = r*sin(u)
plot3d(x,y,z,col = 'blue',alpha = 0.3,type = 'l')
points3d(cos(pi/4)*cos(pi/4),cos(pi/4)*sin(pi/4),sin(pi/4),col = 'red',size = 5)


install.packages("pavo")

##################
###静态四面体
library(pavo)
set.seed(1612217)

specs <- readRDS(system.file("extdata/specsdata.rds",
                             package = "pavo"
))
mspecs <- aggspec(specs, by = 3, FUN = mean)
spp <- gsub("\\.[0-9].*$", "", names(mspecs))[-1]
sppspec <- aggspec(mspecs, by = spp, FUN = mean)
spec.sm <- procspec(sppspec, opt = "smooth", span = 0.2)

data(flowers)

head(flowers[1:4])
vis.flowers <- vismodel(flowers,
                        visual = "bluetit",
                        qcatch = "fi",
                        scale = 10000
)

tetra.flowers <- colspace(vis.flowers, space = "tcs")
head(tetra.flowers)
plot(tetra.flowers,
     pch = 15,
     bg = spec2rgb(flowers),
     #perspective = TRUE,
     range = c(5, 10),
     cex = 0.1
)

##############
###构建四面体###
tetrahedron <-
  rbind(
    c(2*sqrt(2)/3, 0, -1/3),
    c(-sqrt(2)/3, sqrt(2/3), -1/3),
    c(-sqrt(2)/3, -sqrt(2/3), -1/3),
    c(0, 0, 1)
  )
#install.packages("uniformly")
library(uniformly)
set.seed(314)
randomPoints <- runif_in_tetrahedron(
  3, tetrahedron[1, ], tetrahedron[2, ], tetrahedron[3, ], tetrahedron[4, ]
)

plotTetrahedron <- function(tetrahedron, alpha){
  faces <- combn(4L, 3L)
  for(j in 1L:ncol(faces)){
    triangles3d(tetrahedron[faces[, j], ], color = "orange", alpha = alpha)
  }
  edges <- combn(4L, 2L)
  for(j in 1L:ncol(edges)){
    shade3d(
      cylinder3d(tetrahedron[edges[, j], ], sides = 60, radius = 0.03),
      color = "yellow"
    )
  }
  spheres3d(tetrahedron, radius = 0.05, color = "yellow")
}

open3d(windowRect = c(50, 50, 562, 562))
bg3d("slategray")
plotTetrahedron(tetrahedron, alpha = 0.3)
spheres3d(randomPoints, radius = 0.04, color = "black")
####
library(rgl)
data(iris)
head(iris)

x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width
m <- pet.w <- iris$Petal.Width
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
  
  rgl_add_axes <- function(x, y, z, axis.col = "grey",
                           xlab = "", ylab="", zlab="", show.plane = TRUE, 
                           show.bbox = FALSE, bbox.col = c("#333377","black"))
  { 
    
    lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
    xlim <- lim(x); ylim <- lim(y); zlim <- lim(z);
    rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
    rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
    rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
    
    axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                  c(0, 0, zlim[2]))
    rgl.points(axes, color = axis.col, size = 3)
    
    rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
              adj = c(0.5, -0.8), size = 2)
    
    if(show.plane) 
      xlim <- xlim/1.1; zlim <- zlim /1.1
    rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
               z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
    
    if(show.bbox){
      rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
               emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
               xlen = 3, ylen = 3, zlen = 3) 
    }
  }
  
cols <- get_colors(iris$Species, c("#999999", "#E69F00", "#56B4E9"))
rgl_init()
rgl.points(x, y, z, m, r = 0.2, color = cols) 
rgl_add_axes(x, y, z)

##########
update.packages(ask=FALSE,checkBuilt=TRUE,repos="https://cloud.r-project.org")

library(glmnet)
load("CoxExample.RData")

##########
x <- c(1,2,3,4,5,6)
y <- c(0)*length(x)
for (i in (1:length(x)))
  + if(x[i]==7)(y[i] <- 0)else(y[i] <- 1)
print(y)

##########
dat <- head(iris)
rownames(data) <- data[,1]


##########
print("hello,world")
##########
library(tibble)
a <- c(1,2,3,4)
b <- c(2,3,4,5)
c <- c(5,6,7,8)
data=tibble(a,b,c)
data$a
formula=a~b+c
mode=lm(formula,data=data)
summary(mode)
dim(a)

##########
data <- data.frame(A = c(1, 0, 1, 0, 1), B = c(0, 0, 1, 1, 0), C = c(0, 0, 0, 1, 0), D = c(1, 0, 1, 0, 1), E = c(1, 0, 0, 1, 0))
has_disease <- function(row) {
  if (sum(row) > 0) {
    return(1)
  } else {
    return(0)
  }
}
data$has_disease <- apply(data, 1, has_disease)

##########
install.packages("SeuratData")
InstallData("st×Brain")

#########
rm(list = ls())
a <- read.table("CCLE_exp.csv",sep = "\t",check.names = F)

#########
install.packages('ecodist')
install.packages('philentropy')
library(philentropy)
library(ecodist)
library(phyloseq)
iris
library(philentropy)
P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
data <- rbind(P,Q)
distance(data, method = "euclidean")
getDistMethods()

#########
# 生成向量A
A <- seq(1, 100, by = 2)

# 查看A的长度
length(A)

#########
x
for(i in 1:1000)
  {x= c(x, runif(1, min = 0, max = 1))}
plct(density(x))

#########
A <- (1:100)[1:100 %% 2 !=0]
length(A)
(1:100)[1:100 %% 2 !=1]

#########
data <- 1:10
data
round(sd(data), digits = 2)
sd(data)

########
BiocManager::install ("org.Mm.eg.db")

########
library(devtools)
#install_github("jmzeng1314/GEOmirror")
library("GEOmirror")
#install_github("jmzeng1314/idmap1")
#install_github("jmzeng1314/idmap2")
#install_github("jmzeng1314/idmap3")
library(idmap1)
library(idmap2)
library(idmap3)
library("GEOquery")
#install.packages("meta")
library("meta")
#下载数据和注释
eset<-geoChina('GSE131793')
ann_info6244 <-getIDs('GPL6244')
dim(ann_info)
ann_info6244_1 <- getGEO(GEO = 'GPL6244',destdir = ".")
dim()
ann_info6480 <- getGEO(GEO = 'GPL6480',destdir = ".")
dim(ann_info)
#获取表达矩阵
eset
eset<-eset[[1]]
probe_exp <- exprs(eset) 
dim(probe_exp)
data_use <- as.data.frame(probe_exp)
#获取注释信息
ann_info <- as.matrix(ann_info)
class(ann_info)
Meta(ann_info111)$title

?Meta
anno_geoquery<-Table(ann_info)

##########（r语言）为什么这个function单独输入参数的时候可以得出结果，而输入向量的时候会报错呢
rm(list = ls())
Y <- c(sort(rlnorm(10, meanlog = 3, sdlog = 1)))
n <- length(Y)
h <- function(x,y){
  flag <- 0
  if(x>0 && y>0){
    for (i in 2:n) {
      a <- (Y[i-1]-x)/y
      b = (Y[i]-x)/y
      flag <- flag + log(integrate(function(t) t^2 ,a,b)$value)
    }}
  return(flag)
}
h(1,1)
h(2,1)
m <- c(1,2)
h(m,1)
print(m,1)
print(m)


#########
rm(list = ls())
buydates = as.Date(c("2013-07-01","2013-07-09","2013-07-29","2018-08-20"))
sells = as.Date(c("2013-07-02","2013-07-10","2013-07-30","2018-08-21"))

as.Date(ifelse(length(buydates)>length(sells),
               c(sells,as.Date("2023-01-30",origin="1970-01-01")),sells), origin="1970-01-01")
as.Date(sells)


########矩阵内列提取特定字符
rm(list=ls())
# 创建示例数据
data <- data.frame(col = c("abc(def)gh", "ijk(lmn)op", "qrs(tuv)wx"))
# 提取括号内的字符
data$extracted <- gsub(".*\\((.*)\\).*", "\\1", data$col)
# 输出结果
print(data$extracted)
###答案
gene$names <- gsub(".*\\((.*)\\).*", "\\1", gene3$SPOT_ID.1)
print(gene$names)

#######设置断轴
rm(list = ls())
library(ggplot2)
# 创建示例数据
data <- data.frame(category = c("A", "B", "C", "D"), value = c(10, 20, 30, 40))
# 绘制基本的条形图
p <- ggplot(data, aes(x = category, y = value)) + 
  geom_bar(stat = "identity") + 
  labs(x = "Category", y = "Value")
p
# 设置断轴
p + scale_x_continuous(limits = c(1, 4.5), expand = c(0, 0))
# 生成数据
data <- data.frame(date = c("2021-11-18", "2021-12-01", "2021-12-02", "2021-12-03", "2021-12-04", "2021-12-05", "2021-12-06", "2021-12-07", "2021-12-08", "2021-12-09", "2021-12-10", "2021-12-11", "2021-12-12", "2021-12-13", "2021-12-14", "2021-12-15", "2021-12-16", "2021-12-17", "2021-12-18", "2021-12-19", "2021-12-20", "2021-12-21", "2021-12-22", "2021-12-23", "2021-12-24", "2021-12-25", "2021-12-26", "2021-12-27", "2021-12-28", "2021-12-29", "2021-12-30", "2021-12-31"), value = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32))
# 转换日期格式
data$date <- as.Date(data$date)
# 绘制图形
ggplot(data, aes(x = date, y = value)) +
  geom_point() + scale_x_continuous(limits = c(2021-11-19, 2021-11-30), expand = c(0, 0))
# 创建示例数据
data <- data.frame(category = c("A", "B", "C", "D"), value = c(10, 20, 30, 40))
# 绘制基本的柱状图
p <- ggplot(data, aes(x = category, y = value)) + 
  geom_bar(stat = "identity") + 
  labs(x = "Category", y = "Value")
# 使用scale_x_discrete()函数设置离散比例尺
p + scale_x_discrete() + scale_y_continuous()


######批量修改列内变量
# 创建一个包含原始数据的矩阵
data <- matrix(c("A", "B", "A", "C", "C", "D"), nrow = 3, ncol = 2)
data
data <- apply(data, 2, function(x) {
  ifelse(x == "A", "Male", ifelse(x == "B", "Female", ifelse(x == "C", "tumor", "control")))
})
data


#####r语言生成一个包含100个元素的向量^，元素的取值范围为0~ 1000，并且该向量中前40个元素是偶数，
#####第41～80位置的元素为5的倍数，最后20个元素是3的倍数
rm(list = ls())
# 生成前40个偶数
vec <- seq(0, by = 2, length.out = 40)
# 生成第41到80个5的倍数
vec <- c(vec, seq(41, by = 5, length.out = 40))
# 生成最后20个3的倍数
vec <- c(vec, seq(60, by = 3, length.out = 20))
# 如果向量不足100个元素，使用随机数填充
if (length(vec) < 100) {
  vec <- c(vec, sample(0:1000, size = 100 - length(vec), replace = TRUE))
}
# 打印结果
vec


#####孟德尔随机化运行代码后发现结果为阴性。发现漏斗图底部有两个值距离较远，算是离群值吗？怎么去掉离群值呢？
boxplot(x)
# 计算平均值和标准差
mean_x <- mean(x)
sd_x <- sd(x)
# 找出离群值
outliers <- x[x < mean_x - 3 * sd_x | x > mean_x + 3 * sd_x]
# 从数据集中删除离群值
x_clean <- x[!x %in% outliers]


#####怎么安装MendelR包
install.packages("MendelR", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

#####首列导出设置列名
library(openxlsx)
set.seed(123)
df <- data.frame(matrix(runif(25), nrow = 5))
df
df <- cbind(rownames(df),df);names(df)[1] <- 'ID'
df
write.table(df,file = "df.txt",sep = "\t",row.names = F,col.names = T)
write.xlsx(df,file = "df.xlsx", rowNames = F, colNames = T)

#####走迷宫计算最短距离
# 示例数据
grid <- matrix(c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
  1, 2, 3, 3, 3, 3, 3, 3, 2, 1,
  1, 2, 3, 4, 4, 4, 3, 3, 2, 1,
  1, 2, 3, 4, 5, 4, 3, 3, 2, 1,
  1, 2, 3, 4, 4, 4, 3, 3, 2, 1,
  1, 2, 3, 3, 3, 3, 3, 3, 2, 1,
  1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1
), nrow = 10, ncol = 10)

grid[4:6, 5] <- "*"

# 定义起点和终点
start <- c(2, 2)
end <- c(9, 9)

# 初始化距离矩阵和路径矩阵
dist <- matrix(Inf, nrow = nrow(grid), ncol = ncol(grid))
dist[start[1], start[2]] <- 0
path <- matrix(NA, nrow = nrow(grid), ncol = ncol(grid))

# 定义8个方向
directions <- c(-1, 0, 1)

# 迭代更新距离矩阵和路径矩阵
while (TRUE) {
  # 找到当前未处理的最短距离的格子
  unprocessed <- which(!is.na(dist), arr.ind = TRUE)
  if (length(unprocessed) == 0) break
  current <- unprocessed[which.min(dist[unprocessed])]
}
  # 更新相邻的格子的距离和路径
  for (dx in directions) {
    for (dy in directions) {
      # 排除本身和斜向的情况
      if (dx == 0 & dy == 0) next
      if (abs(dx) == abs(dy)) next
      
      x <- current[1] + dx
      y <- current[2] + dy
      # 排除超出边界和无法通过的情况
      if (x < 1 | x > nrow(grid) | y < 1 | y > ncol(grid)) next
      if (grid[x, y] == "*") next
      
      # 计算新的距离
      new_dist
    }
  }

#####paste有什么用
num <- c("a","b","c","d","e","f","g")
data <- c(1:7)
for (i in num) {
 for (dat in data) {
   write.table(dat,file = paste("1.",i,".txt"),sep = "\t")
 }
}

#####
library(dplyr)

# 筛选出重复患者
duplicated_patients <- patient_data %>% 
  group_by(patient_id) %>% 
  filter(n() > 1)

# 选择缺失值最少的行
least_missing <- duplicated_patients %>% 
  arrange(across(everything(), ~sum(!is.na(.)))) %>% 
  slice(1)

# 选择没有缺失值的行
no_missing <- duplicated_patients %>% 
  filter(rowSums(is.na(.)) == 0) %>% 
  slice(1)

# 合并缺失最少的行和没有缺失的行
merged <- bind_rows(least_missing, no_missing)

#####用edit函数
# 创建一个空的矩阵
data_matrix <- matrix(nrow=10, ncol=5)

# 填充数据
data_matrix[,1] <- 1:10                     # id
data_matrix[,2] <- sample(c("男", "女"), 10, replace=TRUE)  # 性别
data_matrix[,3] <- sample(0:100, 10, replace=TRUE)  # 语文成绩
data_matrix[,4] <- sample(0:100, 10, replace=TRUE)  # 数学成绩
data_matrix[,5] <- sample(0:100, 10, replace=TRUE)  # 英语成绩

# 添加列名
colnames(data_matrix) <- c("id", "gender", "chinese", "math", "english")

a <- data.frame(data_matrix)
edit(a)
class(data_matrix)
