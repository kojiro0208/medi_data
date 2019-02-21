#平均3.56のポワソン分布
x <- c(0:9)
prob <- dpois(x ,lambda = 3.56)
plot(x , prob　,type="b" )

#データが離散型の場合は二項分布かポワソン分布を利用する
#データのあてはまりの良さは最尤推定法で確認できる↓


data <-rpois(500 , 3.5)#ラムダ3.5のポワソン分布に従うデータ
hist(data )

#dpois(数,lambda)の確率
logL <- function(m)sum(dpois(data , m , log = TRUE))
lambda <- seq(2,5,0.1)
#尤度のグラフ
plot(lambda , sapply(lambda , logL) , type = "l")

#-----分布の性質------

# ポワソン分布：非負の離散値（カウントデータ）で上限が特に決まっていない
# 二項分布：非負の離散値（カウントデータ）、分散は平均の関数
# 
# 正規分布：連続値で－∞～＋∞の範囲をとる、平均と分散は無関係
# ガンマ分布：連続値で0～＋∞の範囲をとる、分散は平均の関数

#-----分布の性質------

## GLM
d <- read.csv("http://hosho.ees.hokudai.ac.jp/~kubo/stat/iwanamibook/fig/poisson/data3a.csv")
head(d)
#colnames(d) <-  c("種子数" , "個体サイズ" , "肥料の有無")
hist(d$y,freq = FALSE)
lines(density(d$y))
ggplot(data=d,aes(x = y))+
  geom_density()

#ラムダを個体サイズの関数にする
lamda = exp(α　+ β*X)
log(lamda) = α　+ β*X
#↑　リンク関数という(ロジスティックの場合はロジッと関数)

fit <- glm(y ~ x , data = d , family = poisson(link = "log"))
summary(fit)

#最大対数尤度は
logLik(fit)

lamda <- function(x)exp(1.29 + 0.0756*x)

#新たな個体サイズ
xx <- seq(min(d$x) , max(d$x), length = 100)

ggplot(data = d ,aes(x , y ,  color = f ))+
  geom_point()+
  geom_line(aes( xx, lamda(xx)))
predict(fit , newdata = data.frame(x = xx))

#因子型も勝手にダミー変数化してくれる
fit2 <- glm(y ~ f , data = d , family = poisson(link = "log"))
summary(fit2)#fTはTdダミーという意味

logLik(fit2)#さっきより当てはまりが悪い

#全部の変数を使う
fit_all <- glm(y ~ x + f , data = d , family = poisson)
logLik(fit_all)

#正規線形にした場合、種子の数がマイナスになることがある。
fit_all <- glm(y ~ x + f , data = d )

pois <- dpois(c(1:100) , 50)
rpois <- rpois(c(1:10000) , 50)
plot(c(1:100) , pois , type = "l")
hist(rpois , breaks = c(1:100))

###分布に合わせたモデルが必要！

#===4章====#

###____良いモデルとは？____####
fit_null <- glm(y ~ 1 , data = d , family = poisson(link = "log"))
fit <- glm(y ~ x , data = d , family = poisson(link = "log"))
fit2 <- glm(y ~ f , data = d , family = poisson(link = "log"))

fit_all <- glm(y ~ x + f , data = d )
#dpois(数,lambda)の確率

#最小逸脱度（パラメーターをサンプルのサイズ分だけ使った場合）
sum(log(dpois(d$y , d$y)))*-2

#AIC
d$rand_y <- rpois(100 , 8)
hist(d$rand_y)
fit_null <- glm(rand_y~1 , data = d , family = poisson)
fit <- glm(rand_y~x , data = d , family = poisson)
summary(fit_null)
y_pred <- predict(fit_null , newdata = data.frame(d$median_x))

#予測したモデルでリサンプリングした100セットのデータの最大対数尤度を計算する

diff_l <- c()
for( h in c(1:1000)){
  train_y <- rpois(100 , 8)
  fit <- glm(train_y~1 , family = poisson)
  l_list <- c()
  for(i in c(1:200)){
    test_y <- rpois(100 , 8)
    max_l <- sum(dpois(test_y , exp(fit$coefficients) , log = TRUE))
    l_list <- c(l_list , max_l)
  } 
  meanl <- mean(l_list)
  true_value <- sum(dpois(train_y , exp(fit$coefficients) , log = TRUE))
  diff_l <- c(diff_l,true_value-meanl)
  
}



diff_l <- c()

for(i in 1:100){
  x_true <- rnorm(100,10,10)
  y <- 5*x_true+rnorm(100,0,10)
  fit <- lm(y~x_true)
 
  l_vec <- c()
  for(i in 1:100){
    MSE <- fit$residuals^2 %>% mean() %>% sqrt()
    y_sim <- fit$fitted.values+rnorm(100,0,MSE)
    
    max_l <- dnorm(y_sim, fit$fitted.values , MSE ,log=T)%>% sum()
    l_vec <- c(l_vec ,max_l)
  }
  
  meanl <- mean(l_vec)
  true_l <- logLik(fit)
  diff_l <- c(diff_l,true_l-meanl)
}

plot(density(diff_l))
mean(diff_l)




incom <- data(income)
head(incom)

#===5章====#

###____検定____####




fit1 <- glm(y ~1 , data = d , family = poisson)
fit2 <- glm(y ~x , data = d , family = poisson)

fit1$deviance - fit2$deviance

#ランダムデータの生成と逸脱度の差を評価
get.dv <- function(d){
  n.sample = nrow(d)
  lamda = mean(d$y)
  d$y_rand <- rpois(n.sample , lamda)
  
  fit1 <- glm(y_rand ~1 , data = d , family = poisson)
  fit2 <- glm(y_rand ~x , data = d , family = poisson)
  
  return(fit1$deviance - fit2$deviance)
}

#1000回繰り返す
pb <-  c()
for(i in 1:1000){
  pb <- c(pb , get.dv(d))
}

#尤度比検定
summary(pb)
hist(pb , breaks = 100)
sum(pb>4.0)

#95%信頼区間
quantile(pb , 0.95)

#χ二乗分布で近似
anova(fit1 , fit2 , test = "Chisq" )

#===6章====#

###____ロジスティック____####
d <- read.csv("datalog.csv")

ggplot(data = d[d$f == "T",] , aes(x , y , group=f , color =f))+
  geom_point(size = 3)

logistic(4)
z <- seq(-6,6 , 0.1)
plot(z , logistic(z) , type = "l")
library(MASS)
#ステップワイズAIC

fit.xf <- glm(cbind( y , N-y)~x+f , data=d , family =binomial)
summary(fit.xf)
summary(step(fit.xf))

glm(cbind( y , N-y)~x*f , data=d , family =binomial)
