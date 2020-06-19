install.packages("C50")
install.packages("randomForest")
install.packages("kernlab")
install.packages("caret",dependencies = TRUE)

rm(list=ls())

library(MASS)
library(C50)
library(caret)
library(randomForest)
library(tidyverse)
theme_set(theme_classic())

wdbc <- read.csv("wdbc.csv", header = T)
head(wdbc)
features <- c("radius", "texture", "perimeter", "area", "smoothness", "compactness",
              "concavity", "concave_points", "symmetry", "fractal_dimension")
names(wdbc) <- c("id", "diagnosis", paste0(features,"_mean"), paste0(features,"_se"),
                 paste0(features,"_worst"))
wdbc.data <- wdbc[,c(3:32)]
head(wdbc.data)
wdbc$id
row.names(wdbc.data) <- wdbc$id
wdbc_raw <- cbind(wdbc.data, as.factor(as.numeric(wdbc$diagnosis)-1))
colnames(wdbc_raw)[31] <- "diagnosis"
#Note the last part removes ID as a variable, recodes the diagnosis as either 1 or 0 (deletes a few strays!) and
#retains the ID in the rownames.

head(wdbc_raw)
#the last part removes ID as a variable, recodes the diagnosis as either 1 or 0 (deletes a few strays!) and
#retains the ID in the rownames
# Split data into Train and Test dataset
#the data is split into training (75%) and testing/evaluation (25%) sets
ii <- createDataPartition(wdbc_raw[,31], p=3/4, list=F) ## returns indices for train data 
xTrain <- wdbc_raw[ii,]; yTrain <- wdbc_raw[ii]
xTrain
yTrain
xTest <- wdbc_raw[-ii,]; yTest <- wdbc_raw[-ii]
xTest
yTest
# Fit and plot model
yTrain <- as.factor(yTrain)
mdl <- C5.0(x=xTrain, y=yTrain)
plot(mdl)
# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
confusionMatrix(yTestPred, yTest) # predicted/true

ii <- createDataPartition(iris[, 5], p=.7, list=F) ## returns indices for train data
xTrain <- iris[ii, 1:4]; yTrain <- iris[ii, 5]
xTest <- iris[-ii, 1:4]; yTest <- iris[-ii, 5]

# Question 1
classification <- read.csv("Classification.csv")
head(classification)
X1 <- classification[,-2]
head(X1)
X2 <- classification[,-1]
head(X2)
summary(X1)
summary(X2)

#Question 2
# Define function that will split the data into training and test sets
trainTestSplit <- function(df,trainPercent,seed1){
  ## Sample size percent
  smp_size <- floor(trainPercent/100 * nrow(df))
  ## set the seed
  set.seed(seed1)
  train_ind <- sample(seq_len(nrow(df)), size = smp_size)
  train_ind
}

# Split as training and test sets
train_ind <- trainTestSplit(classification,trainPercent=80,seed=123)
train_ind
nrow(classification)
seq_len(nrow(classification))
train.data <- classification[train_ind, ]
train.data
test.data <- classification[-train_ind, ]
test.data

# Question 3
### Linear discriminant analysis.
head(classification)
# Normalize the data
# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

head(classification)
# Fit the model
model <- lda(Group~., data = train.transformed)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy 0.73
mean(predictions$class==test.transformed$Group)
model <- lda(Group~., data = train.transformed)
model
plot(model)
# Make predictions
predictions <- model %>% predict(test.transformed)
names(predictions)
#Inspect the results
# Predicted classes
head(predictions$class, 6)
# Predicted probabilities of class memebership.
head(predictions$posterior, 6) 
# Linear discriminants
head(predictions$x, 3) 
mean(predictions$class==test.transformed$Group)
sum(predictions$posterior[ ,1] >=.5)

### Quadratic discriminant analysis.
# Fit the model
model <- qda(Group~., data = train.transformed)
model
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy 0.79
mean(predictions$class == test.transformed$Group)

### Logistic regression. well done # Lab 7
# Fit a generalized linear logistic model
# Note that the last column (Group) contains the labels
train.data <- classification[train_ind, ]
train.data
test.data <- classification[-train_ind, ]
test.data
fit <- glm(Group~.,family=binomial,data=train.data,control = list(maxit = 50))
# Next we check the performance of the trained model. Here we use the threshold p∗ = 0.5.
# Predict the output from the model
a=predict(fit,newdata=train.data,type="response")
# Set response >0.5 as 1 and <=0.5 as 0
b=ifelse(a>0.5,1,0)
# Compute the confusion matrix for training data
confusionMatrix(as.factor(b),as.factor(train.data$Group)) #Accuracy : 0.7275 

# Finally we use the trained model to make predictions for the test dataset. 
# Again, we use the threshold p∗ = 0.5.
m=predict(fit,newdata=test.data,type="response")
n=ifelse(m>0.5,1,0)
# Compute the confusion matrix for test output
confusionMatrix(as.factor(n),as.factor(test.data$Group))  # Accuracy : 0.73     


### Support vector machines. # Lab 9
# Split test/train
set.seed(103) # for reproducibility
head(classification)
classification$Group=factor(classification$Group,levels = c(0,1))
ii <- createDataPartition(classification[,3],p=.8,list = F)
xTrain <- classification[ii, 1:2]; yTrain <- classification[ii, 3]
xTest <- classification[-ii, 1:2]; yTest <- classification[-ii, 3]
dim(xTrain)
dim(xTest)
length(yTest)
head(yTest)
library(kernlab)
library(e1071)
library(ISLR)
library(RColorBrewer)


mdl <- train(x=xTrain,y=yTrain, method='svmLinear')
print(mdl)
# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
yTestPred
#yTestPred <- mdl %>% predict(xTest)
table(yTestPred, yTest)
confusionMatrix(yTestPred, yTest) # predicted/true Accuracy : 0.7085   

mdl <- train(x=xTrain,y=yTrain, method = "svmLinear",
             trControl = trainControl("cv", number = 5),
             tuneGrid = expand.grid(C = seq(0, 2, length = 20)))

# Plot model accuracy vs different values of Cost
plot(mdl)
# Print the best tuning parameter C that maximises model accuracy
mdl$bestTune
# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
confusionMatrix(yTestPred, yTest) # predicted/true Accuracy : 0.7085 

mdl <- train(x=xTrain,y=yTrain, method='svmPoly')
print(mdl)
mdl$bestTune
# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
confusionMatrix(yTestPred, yTest) # predicted/true  

### K-nearest neighbour regression.
# Set training options
# Repeat 5-fold cross-validation, ten times
opts <- trainControl(method='repeatedcv', number=5, repeats=10, p=0.7)
# Find optimal k (model)
set.seed(1040) # for reproducibility
mdl <- train(x=xTrain, y=yTrain, # training data
             method='knn', # machine learning model
             trControl=opts, # training options
             tuneGrid=data.frame(k=seq(2, 15))) # range of k's to try
print(mdl)

# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
length(yTestPred)
yTestPred <- as.factor(yTestPred)
yTest <- as.factor(yTest)
confusionMatrix(yTestPred, yTest) # predicted/true  Accuracy : 0.7688



