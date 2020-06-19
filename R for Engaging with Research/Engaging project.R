install.packages("quantmod")
installed.packages("TTR")
installed.packages("zoo")
installed.packages("xts")
install.packages("PerformanceAnalytics")
install.packages("scales")

library(quantmod)
library(TTR)
library(PerformanceAnalytics)
library(ggplot2)
library(scales)
# Get data 
options(stringsAsFactors = FALSE)
symbols <- c("^GSPC","^N225","^HSI","^STI","000001.SS")
suppressWarnings(getSymbols(symbols,src = "yahoo",from="2012-01-01"))

df <- merge(GSPC$GSPC.Adjusted,HSI$HSI.Adjusted,N225$N225.Adjusted,STI$STI.Adjusted,
            `000001.SS`$`000001.SS.Adjusted`)
names(df) <- c("GSPC","HSI","N225","STI","SSE")

head(df)

# Get stock data
setSymbolLookup(SS=list(name='000001.ss',src='yahoo',from='2015-01-01',to='2020-01-01'))
getSymbols("SS")
chartSeries(SS)
setSymbolLookup(WK=list(name='000002.sz',src='yahoo'))
getSymbols("WK")
chartSeries(WK)

# Get 000001.SZ data
options(stringsAsFactors = FALSE)
symbols <- c("000001.SZ")
suppressWarnings(getSymbols(symbols,src = "yahoo",from="2018-01-01"))
df <- merge(`000001.SZ`$`000001.SZ.Open`,`000001.SZ`$`000001.SZ.High`,`000001.SZ`$`000001.SZ.Low`,`000001.SZ`$`000001.SZ.Close`,`000001.SZ`$`000001.SZ.Volume`)
names(df) <- c("Open","High","Low","Colse","Volume")
head(df)
class(df)

# Get stock data
getSymbols(c('IBM', 'BABA', '^GSPC','NTES'), from = '2014-09-19')
names(BABA)
head(BABA)

names(BABA) = c('open', 'high', 'low', 'close', 'volume', 'adjusted')
SP500 = GSPC
names(SP500) = c('open', 'high', 'low', 'close', 'volume', 'adjusted')
names(IBM) = c('open', 'high', 'low', 'close', 'volume', 'adjusted')
names(NTES) = c('open', 'high', 'low', 'close', 'volume', 'adjusted')

BABAdf = as.data.frame(BABA)
str(BABAdf)
head(BABAdf) 
BABAdf['2014-09-19',]

SP500df = as.data.frame(SP500)
IBMdf = as.data.frame(IBM)
NTESdf = as.data.frame(NTES)

stock = list(SP500 = SP500df, IBM  = IBMdf, NTES = NTESdf, BABA = BABAdf)
str(stock)
head(stock$SP500)
head(stock[['SP500']])

dateArea <- function(sDate = Sys.Date()-365, 
                     eDate = Sys.Date(), before = 0){
  if(class(sDate)=='character') sDate = as.Date(sDate)
  if(class(eDate)=='character') eDate = as.Date(eDate)
  return(paste(sDate - before, eDate, sep= '/'))
}

ma <- function(cdata, mas = c(5, 20, 60)){
  if(nrow(cdata) <= 60) return(NULL)
  ldata <- cdata
  for(m in mas){
    ldata <- merge(ldata, SMA(cdata, m))
  }
  names(ldata) <- c('value', paste('ma', mas, sep = ''))
  return(ldata)
}

title = 'BABA'
sDate <- as.Date('2015-01-01')
eDate <- as.Date('2015-3-20')
cdata <- BABA[dateArea(sDate, eDate, 360)]$close
ldata = NULL
ldata <- ma(cdata, c(5, 20, 60))
tail(ldata)
#get K line
plot(ldata)
getMaSd <- function(ldata, mas = 20, sDate, eDate){
  if(is.null(ldata) || nrow(ldata) <= max(mas)) return(NULL)
  col <- paste('ma', mas,sep = '')
  ldata <- ldata[,c('value', col)]
  ldata$dif <- ldata[,col] - ldata$value
  ldata$sd <- runSD(ldata[,'dif'], mas)
  ldata$rate<- round(ldata$dif/ldata$sd, 2)
  ldata[dateArea(sDate, eDate)]
}

ldata5 <- getMaSd(ldata, 5, sDate, eDate)
ldata20 <- getMaSd(ldata, 20, sDate, eDate)
ldata60 <- getMaSd(ldata, 60, sDate, eDate)
head(ldata20)
plot(ldata5[,c('dif', 'sd', 'rate')])
plot(ldata20[,c('dif', 'sd', 'rate')])
plot(ldata60[,c('dif', 'sd', 'rate')])

#find the poibnt that has sd twice greater than the average sd
#dir = 1, buy; dir = 2, short
buypoint <- function(ldata, x = 2, dir = 2){
  idx <- which(ldata$rate >x)
  if(dir ==2){
    idx <- c(idx, which(ldata$rate<x*-1))
  }
  return(ldata[idx,])
}
buypoint5 <- buypoint(ldata5)
str(buypoint5)

plot(ldata5[,c('dif', 'sd', 'rate')])
points(buypoint5[,'rate'])
buypoint20 <- buypoint(ldata20, dir =1)
plot(ldata20[,c('dif', 'sd', 'rate')])
points(buypoint20[,'rate'])
points(buypoint20[,'sd'], col = 'red')

#calculate the selling points
sellpoint <- function(ldata, buydata){
  buy = buydata[which(buydata$dif >0),]
  
  aidx <- index(ldata[which(ldata$dif <= 0),])
  sellidx <- sapply(index(buy), function(ele){
    head(which(aidx >ele), 1)
  })
  ldata[aidx[unique(unlist(sellidx))]]
}

sellpoint5 <- sellpoint(ldata5, buypoint5)
head(sellpoint5)
plot(ldata5[,])
points(buypoint5[,'value'])
points(sellpoint5[,'value'], col = 'red')

bsdata <- merge(buypoint5$value, sellpoint5$value)
head(bsdata)
names(bsdata) = c('buy', 'sell')
plot(ldata5[,'value'])
points(bsdata[,'buy'], col = 'red')
points(bsdata[,'sell'], col='green')

signal <- function(buy, sell){
  selldf <- data.frame(sell, op = as.character(rep('S', nrow(sell))))
  buydf <- data.frame(buy, op = as.character(rep('B', nrow(buy))))
  sdata <- rbind(buydf, selldf)
  sdata[order(as.Date(row.names(sdata))),]
}

sdata <- signal(buypoint5, sellpoint5)
sdata
#the default setting for this trade is
# capital 100000, every buy action costs 10000
trade <- function(sdata, capital = 100000, fixMoney = 10000){
  amount = 0
  cash = capital
  ticks <- data.frame()
  
  for(i in 1:nrow(sdata)){
    row <- sdata[i,]
    if(row$op == 'B'){
      if(cash < fixMoney){
        print(paste(row.names(row), 'Not enough cash'))
        next
      }
      amount0<- floor(fixMoney/row$value)
      amount <- amount + amount0
      cash <- cash - amount0*row$value
    }
    if(row$op =='S'){
      cash <- cash+amount*row$value
      amount = 0
    }
    
    row$cash <- round(cash,2)
    row$amount <- amount
    row$asset <- round(cash+amount*row$value, 2)
    ticks <- rbind(ticks, row)
  }
  
  ticks$diff <- c(0, round(diff(ticks$asset), 2))
  rise <- ticks[intersect(which(ticks$diff >0), which(ticks$op =='S')),]
  fall <- ticks[intersect(which(ticks$diff <0), which(ticks$op =='S')),]
  return(list(ticks = ticks, rise = rise, fall = fall))
}

return5 <- trade(sdata)

str(return5)

head(return5$ticks)
head(return5$rise)
head(return5$fall)

class(return5$rise)

tick5 <- return5$ticks
str(tick5)
head(tick5)
index(tick5)
tick5$date = as.Date(row.names(tick5))

head(tick5)

tick50 <- zoo(tick5[,'value'], as.Date(tick5$date))
class(tick50)
str(tick50)

tick50 <- as.xts(tick50)
names(tick50) = 'value'
plot(tick50)

names(tick50)

head(tick50$'value')

tick50$asset = tick5[,'asset']
tick50$op = tick5[,'op']
head(tick50)

par(mfrow = c(2,1))
plot(ldata5[,'value'], main = 'stock price of BABA')
points(tick50[which(tick50$'op' == 1), 'value'], col ='red')
points(tick50[which(tick50$'op' == 2), 'value'], col ='green')
plot(tick50[,'asset'], main = 'the increment in assets')


#try all the functions with other stocks
#cdata <- BABA[dateArea(sDate, eDate, 360)]$close; ldata = NULL; ldata <- ma(cdata, c(5, 20, 60))
quick <- function(stock.data, sDate, eDate){
  cdata <- stock.data[dateArea(sDate, eDate, 360)]$close
  ldata <- ma(cdata, c(20))
  ldata <- getMaSd(ldata, 20, sDate, eDate)
  buydata <- buypoint(ldata, 2, 1)
  selldata <- sellpoint(ldata, buydata)
  sdata <- signal(buydata, selldata)
  return(trade(sdata))
}

babaop <- quick(BABA, sDate, eDate)
ntesop <- quick(NTES, '2017-12-01', '2018-02-20')

ntesop$ticks
class(ntesop$ticks)
ntes.tick = zoo(ntesop$ticks[,c('value', 'asset')], as.Date(row.names(ntesop$ticks)))
head(ntes.tick)
ntes.tick$'op' = ntesop$ticks[,'op']
head(ntes.tick)
ntes.tick <- as.xts(ntes.tick)

cdata <- NULL
cdata <- NTES[dateArea('2019-12-01', '2020-05-20', 360),]$close
ldata <- ma(cdata, c(20))
ldata <- getMaSd(ldata, 20, '2019-12-01', '2020-05-20')
head(ldata)
par(mfrow = c(2,1))
plot(ldata[,c('value', 'ma20')], main = 'stock price of NTES')
points(ldata[index(ntes.tick[which(ntes.tick$'op'==1),]),'value' ], col ='red')
points(ldata[index(ntes.tick[which(ntes.tick$'op'==2),]),'value'], col ='green')

plot(ntes.tick[,'asset'], main = 'the increment in assets')

