classification <- read.csv("ClassificationTrue.csv")
head(classification)
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
# Model accuracy 0.73 -> 0.795
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
confusionMatrix(yTestPred, yTest) # predicted/true   

mdl <- train(x=xTrain,y=yTrain, method = "svmLinear",
             trControl = trainControl("cv", number = 5),
             tuneGrid = expand.grid(C = seq(0, 2, length = 20)))

# Plot model accuracy vs different values of Cost
plot(mdl)
# Print the best tuning parameter C that maximises model accuracy
mdl$bestTune
# Test model on testing data
yTestPred <- predict(mdl, newdata=xTest)
confusionMatrix(yTestPred, yTest) # predicted/true 

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
confusionMatrix(yTestPred, yTest) # predicted/true  



