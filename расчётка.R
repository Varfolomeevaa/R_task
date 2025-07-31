#1. предобработка
data = RUS_INVFC_M #инвестиции в основной капитал
data = setNames(data,c('amount','year'))
year = as.numeric(sub("\\..*", "", data$year))
month = as.numeric(sub(".*\\.", "", data$year))
month = ifelse(nchar(sub(".*\\.", "", data$year)) == 1,month * 10, month)
data$year = as.Date(paste(year,month,'01',sep="-"))
plot(x=data$year,y=data$amount,type='l')
plot(x=data$year,y=log(data$amount),type='l')
data$amount = log(data$amount)

#2.а обработка тренда
t = 1:(NROW(data)-10)
fit1=lm(data$amount[1:230]~t) 
summary(fit1)
plot(data$amount,type='l')
lines(fit1$fitted.values,col='red') 
plot(fit1$residuals,type='l') 
AIC(fit1)

t2 = t^2
fit2=lm(data$amount[1:230]~t + t2)
summary(fit2)
plot(data$amount,type='l')
lines(fit2$fitted.values,col='red')
plot(fit2$residuals,type='l')
AIC(fit2)

t3 = sqrt(t)
fit3=lm(data$amount[1:230]~t+t3)
summary(fit3)
plot(data$amount,type='l')
lines(fit3$fitted.values,col='red')
plot(fit3$residuals,type='l')
AIC(fit3)

t4 = t^(1/3)
fit4=lm(data$amount[1:230]~t+t4)
summary(fit4)
plot(data$amount,type='l')
lines(fit4$fitted.values,col='red')
plot(fit4$residuals,type='l')
AIC(fit4)

t5 = t ^ 3
fit5 = lm(data$amount[1:230]~t+t2+t5)
summary(fit5)
plot(data$amount,type='l')
lines(fit5$fitted.values,col='red')
plot(fit5$residuals,type='l')
AIC(fit5)

#2.б остановимся на 4 модели,анализируем сезонность
s = array(1:12,dim=230) 
fit=lm(data$amount[1:230]~t + t4+ factor(s)) 
summary(fit)
anova(fit)
plot(data$amount,type='l')
lines(fit$fitted.values,col='red') 
plot(fit$residuals,type='l')
AIC(fit)

#2.в интервальный прогноз
p3 = predict(fit,data.frame(t=231:240, t4 = (231:240) ^ (1/3),s =array(1:12,dim=10)),interval='prediction')
p3
plot(data$amount,type='l', xlim = c(1,240), ylim=range(c(data$amount,p3)))
lines(fit$fitted.values,col='red')
lines(x = 231:240, y = p3[,1], col = "blue")
lines(x = 231:240, y = p3[,2],col='green') #доверительные интервалы
lines(x = 231:240, y = p3[,3],col='green')

#3. проверка на стационарность
acf(fit$residuals) # нестационарность
r=cor(fit$residuals,t,method='spearman')
r*sqrt((230-2)/(1-r^2))
# t238 = 1.96 расчётное значение меньше статистики, поэтому Н0 принимаем, нет тренда, процесс случайный
# тест дики-фулера
library(forecast)
ndiffs(fit$residuals,test='adf') #стационарность

#4.а полиномиальное сглаживание
pma=function(data_,p=1,m=10){
  result=data_
  t=(-m):m
  for(i in (m+1):(NROW(data_)-m)){
    fit_ = switch(paste(p), #if. paste - текст
                 '1'=lm(data_[(i-m):(i+m)]~t),
                 '2'=lm(data_[(i-m):(i+m)]~t+I(t^2)),
                 '3'=lm(data_[(i-m):(i+m)]~t+I(t^2)+I(t^3)))
    result[i] = fit_$coefficients[1]
    if(i==(m+1)){
      result[1:i]=predict(fit_, newdata=data.frame(t=(-m):0))
    }
    if(i==(NROW(data_)-m)){
      result[i:NROW(data_)]=predict(fit_, newdata=data.frame(t=0:m))
    }
  }
  
  return(result)
}
result1 = pma(data$amount,m=10,p=1)  
plot(data$amount,type='l')
lines(result1,col='red')
result2 = pma(data$amount,m=20,p=3)  
lines(result2,col='blue') #описывает динамику лучше

#4.б экспоненциальное сглаживание
ema=function(data_, alpha=0.5){
  result=data_
  for (i in 2:NROW(data_)){
    result[i]=alpha*data_[i]+(1-alpha)*result[i-1]  
  }
  return(result)
}
result1=ema(data$amount,alpha=0.2)
plot(data$amount,type='l')
lines(result1,col='red')
result2=ema(data$amount,alpha=0.4)
lines(result2,col='blue')
result3=ema(data$amount,alpha=0.7)
lines(result3,col='green') #описывает динамику лучше

#5. Спектральный и гармонический анализ

#5.а разложение в ряд Фурье
res = fit4$residuals 
alpha = c(0:115) 
beta = c(0:115) 
beta[1] = 0
beta[116] = 0
alh=function(res,alpha){
  for (j in 0:115){
    if (j == 0){
      alpha[1] = mean(res)
    }else if (j == 115){
      result = 0
      for (t in 1:230){
        result = result + res[t] * (-1)^t
      }
      alpha[116] = (1 / 230) * result
    } else {
      result = 0
      for (t in 1:230){
        result = result + res[t] * cos(2*pi*j*t/230)
      }
      alpha[j+1] = (2/230) * result
    }
  }
  return(alpha)
}
alpha = alh(res,alpha)
alpha
bet=function(res,beta){
  for (j in 1:114){
    result = 0
    for (t in 1:230){
      result = result + res[t] * sin(2*pi*j*t/230)
    }
    beta[j+1] = (2/230) * result
  }
  return(beta)
}
beta=bet(res,beta)

ph = c(1:230)
pht=function(ph,alpha,beta){
  for (t in 1:230){
    result=0
    for (j in 0:115){
      result = result + alpha[j+1] * cos(2*pi*j*t/230) + beta[j+1]*sin(2*pi*j*t/230)
    }
    ph[t] = result
  }
  return(ph)
}
ph=pht(ph,alpha,beta)

#5.б и в выборочный спектр
s1=spectrum(res)
f1=s1$freq[s1$spec==max(s1$spec)]
f1
1/f1 
fit1=lm(data$amount[1:230]~t+t4+
          cos(2*pi*f1*t)+sin(2*pi*f1*t))
summary(fit1)
AIC(fit1)
plot(data$amount,type='l')
lines(fit1$fitted.values,col='red') 

s2=spectrum(fit1$residuals)
f2=s2$freq[s2$spec==max(s2$spec)]
f2
1/f2
fit2=lm(data$amount[1:230]~t+t4+
          cos(2*pi*f1*t)+sin(2*pi*f1*t)+
          cos(2*pi*f2*t)+sin(2*pi*f2*t))
summary(fit2)
AIC(fit2)
plot(data$amount,type='l')
lines(fit2$fitted.values,col='red')

s3=spectrum(fit2$residuals)
f3=s3$freq[s3$spec==max(s3$spec)]
f3
1/f3
fit3=lm(data$amount[1:230]~t+t4+
          cos(2*pi*f1*t)+sin(2*pi*f1*t)+
          cos(2*pi*f2*t)+sin(2*pi*f2*t)+
          cos(2*pi*f3*t)+sin(2*pi*f3*t))
summary(fit3)
AIC(fit3)
plot(data$amount,type='l')
lines(fit3$fitted.values,col='red')

s4 = spectrum(fit3$residuals)
f4=s4$freq[s4$spec==max(s4$spec)]
f4
1/f4
fit4=lm(data$amount[1:230]~t+t4+
          cos(2*pi*f1*t)+sin(2*pi*f1*t)+
          cos(2*pi*f2*t)+sin(2*pi*f2*t)+
          cos(2*pi*f3*t)+sin(2*pi*f3*t)+
          cos(2*pi*f4*t)+sin(2*pi*f4*t))
summary(fit4)
AIC(fit4)
plot(data$amount,type='l')
lines(fit4$fitted.values,col='red')

s5 = spectrum(fit4$residuals)
f5=s5$freq[s5$spec==max(s5$spec)]
f5
1/f5


acf(fit5$residuals) # стационарность?
r=cor(fit5$residuals,t,method='spearman')
r*sqrt((230-2)/(1-r^2))
# t238 = 1.96 расчётное значение меньше статистики, поэтому Н0 принимаем, нет тренда, процесс случайный
# тест дики-фулера
ndiffs(fit5$residuals,test='adf') #стационарность

#6.а,б,в,г arma
data_arma=fit5$residuals[1:220]
model1 = arima(data_arma,order=c(1,0,1))
summary(model1)

model2 = arima(data_arma,order=c(2,0,1))
summary(model2)

model_auto=auto.arima(data_arma)
summary(model_auto)

plot(forecast(model_auto,10))
lines(x = 221:230, fit5$residuals[221:230], col='red')

#7. тренд+arma+garch
library(fGarch)
gfit1=garchFit(formula = ~ garch(1,1),mode_auto$residuals) 
summary(gfit1)
gfit2=garchFit(formula = ~ garch(1,3),model_auto$residuals) 
summary(gfit2)
p=predict(gfit1,10)
plot(p$standardDeviation,type='l')

#8.а коинтеграция
data1 = BEL_Food
data2 = BGR_Food

data1 = setNames(data1,c('amount','year'))
year = as.numeric(sub("\\..*", "", data1$year))
month = as.numeric(sub(".*\\.", "", data1$year))
month = ifelse(nchar(sub(".*\\.", "", data1$year)) == 1,month * 10, month)
data1$year = as.Date(paste(year,month,'01',sep="-"))
data1

data2 = setNames(data2,c('amount','year'))
year = as.numeric(sub("\\..*", "", data2$year))
month = as.numeric(sub(".*\\.", "", data2$year))
month = ifelse(nchar(sub(".*\\.", "", data2$year)) == 1,month * 10, month)
data2$year = as.Date(paste(year,month,'01',sep="-"))
data2

plot(x=data1$year,y=data1$amount,type='l',ylim=range(c(data1$amount,data2$amount)))
lines(x=data2$year,y=data2$amount,type='l',col='blue')

#8. б определение порядка интегрирования
ndiffs(data1$amount,test='adf')
ndiffs(data2$amount[1:167],test='adf')

#8.в определение коинтеграции
fit_coint=lm(data1$amount~data2$amount[1:167])
summary(fit_coint)
ndiffs(fit_coint$residuals,test='adf')
#некоинтегрированы

#9. а деревья решений
library(randomForest)
data_rf = bank_transactions[1:100000,]
data_rf = na.omit(data_rf) #удаление строк с пропусками
data_rf$TransactionID = NULL
data_rf$CustomerID = NULL
data_rf$CustomerDOB = as.Date(sub(" UTC", "", data_rf$CustomerDOB))
library(eeptools)
data_rf = data_rf[data_rf$CustomerDOB <= Sys.Date(), ]
data_rf
data_rf$CustomerAge = floor(age_calc(
  data_rf$CustomerDOB, 
  enddate = Sys.Date(), 
  units = "years"
))
data_rf$CustomerDOB = NULL
data_rf$CustGender = ifelse(data_rf$CustGender == 'M',1,0)
length(unique(data_rf$CustLocation))
data_rf$TransactionDate = NULL
names(data_rf)[names(data_rf) == "TransactionAmount (INR)"] = "TransactionAmount"
data_rf$CustLocation = NULL
train=data_rf[1:80000,]
test = data_rf[80001:NROW(data_rf),]

ndiffs(train$TransactionAmount,test='adf') #стацион
rf1=randomForest(TransactionAmount~.,data=train,ntree=10)
rfpr1 = predict(rf1,test)
rf2=randomForest(TransactionAmount~.,data=train,ntree=15)
rfpr2 = predict(rf2,test)
rf3=randomForest(TransactionAmount~.,data=train,ntree=20)
rfpr3 = predict(rf3,test)
sum((test$TransactionAmount-rfpr1)^2)
sum((test$TransactionAmount-rfpr2)^2)
sum((test$TransactionAmount-rfpr3)^2)
#выбираем 15 деревьев
importance(rf2) 

fit = lm(TransactionAmount~.,data=train)
lmpr=predict(fit,test)
sum((test$TransactionAmount-lmpr)^2)

#9. б Кластеризация
data_cl = data_rf[1:1000,c('CustAccountBalance','TransactionAmount')]
plot(data_cl$TransactionAmount,data_cl$CustAccountBalance)
#стандартизуем
data_cl$CustAccountBalance = data_cl$CustAccountBalance /sd(data_cl$CustAccountBalance )
data_cl$TransactionAmount= data_cl$TransactionAmount/sd(data_cl$TransactionAmount)

data_cl = data_cl[data_cl$TransactionAmount<4 & data_cl$CustAccountBalance<10,]
plot(data_cl$TransactionAmount,data_cl$CustAccountBalance)

centr = data.frame(
  Balance = c(
    min(data_cl$CustAccountBalance),
    median(data_cl$CustAccountBalance),
    max(data_cl$CustAccountBalance)
  ),
  Amount = c(
    median(data_cl$TransactionAmount),
    min(data_cl$TransactionAmount),
    max(data_cl$TransactionAmount)
  )
)
points(centr$Amount,centr$Balance,col='red')   
centr
data_cl$kl=NA

for (i in 1:NROW(data_cl)){
  x = (data_cl$TransactionAmount[i]-centr$Amount)^2+(data_cl$CustAccountBalance[i]-centr$Balance)^2
  data_cl$kl[i] = (1:3)[x==min(x)] # усл
}                                                                                              
data_cl
plot(data_cl$TransactionAmount,data_cl$CustAccountBalance)
for (i in 1:3){
  points(data_cl$TransactionAmount[data_cl$kl==i],data_cl$CustAccountBalance[data_cl$kl==i],col=rainbow(3)[i]) #генерит 3 цвета
}
centr = aggregate(list(Balance=data_cl$CustAccountBalance, 
                       Amount=data_cl$TransactionAmount), 
                   by=list(kl=data_cl$kl), 
                   mean)
centr

diff = 1
while(diff>0){
  temp = data_cl$kl
  data_cl$kl=NA
  for (i in 1:NROW(data_cl)){
    x = (data_cl$TransactionAmount[i]-centr$Amount)^2+(data_cl$CustAccountBalance[i]-centr$Balance)^2
    data_cl$kl[i] = (1:3)[x==min(x)]
  }                                                                                              
  for (i in 1:3){
    points(data_cl$TransactionAmount[data_cl$kl==i],data_cl$CustAccountBalance[data_cl$kl==i],col=rainbow(3)[i]) #генерит 3 цвета
  }
  centr = aggregate(list(Balance=data_cl$CustAccountBalance, 
                         Amount=data_cl$TransactionAmount), 
                    by=list(kl=data_cl$kl), 
                    mean)
  diff = sum((temp-data_cl$kl)^2)
}
plot(data_cl$TransactionAmount,data_cl$CustAccountBalance)
for (i in 1:3){
  points(data_cl$TransactionAmount[data_cl$kl==i],data_cl$CustAccountBalance[data_cl$kl==i],col=rainbow(3)[i]) #генерит 3 цвета
}

km = kmeans(data_cl,3)

plot(x = data_cl$TransactionAmount, 
     y = data_cl$CustAccountBalance,
     col = km$cluster)

