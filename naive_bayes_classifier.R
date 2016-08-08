# load required packages
library(dplyr)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(e1071)

#seet ggplot2 theme
theme = theme_set(theme_minimal())
theme = theme_update(legend.position="top")

#Richard Shanahan  
#https://github.com/rjshanahan  
#6 August 2016

###### INFS 5094: Assignment 1 Breast Cancer Wisconsin

###### 1. read in file ###### 
breast_cancer <- read.csv('breast-cancer-wisconsin.csv',
                              header=T,
                              sep=",",
                              quote='"',
                              # colClasses=c(
                              #   'character',# id number
                              #   'numeric',   # Clump_Thickness
                              #   'numeric',   # Uniformity_Cell_Size
                              #   'numeric',   # Uniformity_Cell_Shape
                              #   'numeric',   # Marginal_Adhesion
                              #   'numeric',   # S_E_Cell_Size
                              #   'numeric',   # Bare_Nuclei
                              #   'numeric',   # Bland_Chromatin
                              #   'numeric',   # Normal_Nucleoli
                              #   'numeric',  # Mitoses
                              #   'factor'     # class
                              # ),
                              strip.white=T,
                              stringsAsFactors=F,
                              fill=T)

#convert to numeric
breast_cancer$Bare_Nuclei <- as.numeric(breast_cancer$Bare_Nuclei)

#inspect
str(breast_cancer)
describe(breast_cancer)

#check for duplicate records based on ID
nrow(breast_cancer) - length(unique(breast_cancer$id))

#check if there are any missing values
colSums(is.na(breast_cancer)) 

breast_cancer <- breast_cancer %>% filter(!is.na(Bare_Nuclei))

colSums(is.na(breast_cancer)) 

dim(breast_cancer)

#display duplicate id's
dupeid <- breast_cancer$id[duplicated(breast_cancer$id)]
dupeid

dupes <- subset(breast_cancer, duplicated(breast_cancer$id) == T)
dupes <- dupes[order(dupes$id),]

#check for duplicate records based on all equal records
dupeall <- breast_cancer$id[duplicated(breast_cancer)]

dupealldf <- subset(breast_cancer, duplicated(breast_cancer) == T)
dupealldf <- dupealldf[order(dupeall),]

#recode 'class' as character
breast_cancer$class_char <- ifelse(breast_cancer$class == 2,
                                   "benign",
                                   "malignant")



###### 2. visualise distributions ###### 

#class distribution
ggplot(data = breast_cancer, 
       aes(x=class_char,
           fill=class_char)) + 
  geom_bar() +
  ggtitle("UCI Breast Cancer 'Class' Distribution")


#reshape dataset for boxplot representation
breast_cancer.m <- melt(select(breast_cancer, -class, -class_char),
                            id.var="id")

#distribution of other variables
ggplot(data = breast_cancer.m, 
       aes(x=value),
       fill=class) + 
  #geom_histogram(aes(fill=factor(variable)),
  #               binwidth = 1) +
  geom_bar(aes(fill=factor(variable))) +
  facet_wrap("variable") +
  ggtitle("UCI Breast Cancer Dataset Distributions (bin=1)")


#create boxplots of each variable in same graphic
ggplot(data = breast_cancer.m, 
        aes(x=variable, y=value)) + 
        geom_boxplot(aes(fill=variable)) + 
        xlab("Breast Cancer Attribute") + 
        ylab("Value") + 
        ggtitle("UCI Breast Cancer Dataset Boxplots") +
        guides(fill=guide_legend(title="Attribute Legend"))


###### 3. assess correlations ###### 

breast_cancer.cor <- round(cor(breast_cancer[2:10], 
    use = "complete.obs",
    y=NULL,
    method = "pearson"), 2)

breast_cancer.cor


#alternative scatterplot
pairs(breast_cancer[2:10],
      main="Scatterplot of Breast Cancer Numeric",
      pch = 20,
      col="goldenrod")





###### 4. build Naive Bayes classifier - non-dsicretised ###### 

#remove variables not required
breast_cancer_nb <- breast_cancer %>% select(-id, -class)

#split dataset into TRAIN and TEST
breast_cancer_train <- sample_frac(breast_cancer_nb, 0.8)
myrows <- as.numeric(rownames(breast_cancer_train)) 
breast_cancer_test <- breast_cancer_nb[-myrows,]


#inspect split datasets
str(breast_cancer_train)
str(breast_cancer_test)

#summary statistics by class
describe(filter(breast_cancer_train, class_char == 'benign'))
describe(filter(breast_cancer_train, class_char == 'malignant'))


#train naiveBayes model using e1071 package
nb_numeric <- naiveBayes(as.factor(class_char)~.,
                  data=breast_cancer_train,
                  laplace = 0)

#view model summary
nb_numeric
summary(nb_numeric)

#predictions
#prediction function
myPredictor <-function(nb_object, newdata) {
  predict(nb_object,
          newdata=newdata[,-10])  
}

#prediction object
myPredict <- myPredictor(nb_numeric,
                         breast_cancer_test)


#confusion matrix
table(prediction=myPredict,
      true=breast_cancer_test$class_char)

#accuracy
mean(myPredict == breast_cancer_test$class_char)



#10 fold cross-validation
breast_cancer_xval <- tune(naiveBayes,
                           as.factor(class_char)~.,
                           data = breast_cancer_train,
                           predict.func = myPredictor
                           )

#review Cross Validation output
breast_cancer_xval$best.parameters
breast_cancer_xval$best.performance
breast_cancer_xval$performances
breast_cancer_xval$best.model


#use cross validated best model
nb_numeric_xv <- breast_cancer_xval$best.model

#run against test dataset
#prediction object
myPredict_xv <- myPredictor(nb_numeric_xv,
                         breast_cancer_test)


#confusion matrix
table(prediction=myPredict_xv,
      true=breast_cancer_test$class_char)

#accuracy
mean(myPredict_xv == breast_cancer_test$class_char)




###### 5. build Naive Bayes classifier - dsicretised ###### 

#function to dsicretise data into three bins - equal frequency
cutter <- function(myCol, quants) {
  require(Hmisc)
  
  myCol <- cut2(myCol,
                g = quants)
  
}

breast_cancer_nb <- breast_cancer %>% select(-id, -class)

#discretise each variable into 3 equal frequency bins
breast_cancer_nb$Clump_Thickness <- cutter(breast_cancer_nb$Clump_Thickness, 3)
breast_cancer_nb$Uniformity_Cell_Size <- cutter(breast_cancer_nb$Uniformity_Cell_Size, 3)
breast_cancer_nb$Uniformity_Cell_Shape <- cutter(breast_cancer_nb$Uniformity_Cell_Shape, 3)
breast_cancer_nb$Marginal_Adhesion <- cutter(breast_cancer_nb$Marginal_Adhesion, 3)
breast_cancer_nb$S_E_Cell_Size <- cutter(breast_cancer_nb$S_E_Cell_Size, 3)
breast_cancer_nb$Bare_Nuclei <- cutter(breast_cancer_nb$Bare_Nuclei, 3)
breast_cancer_nb$Bland_Chromatin <- cutter(breast_cancer_nb$Bland_Chromatin, 3)
breast_cancer_nb$Normal_Nucleoli<- cutter(breast_cancer_nb$Normal_Nucleoli, 3)
breast_cancer_nb$Mitoses <- cutter(breast_cancer_nb$Mitoses, 3)

#inspect after discretisation
str(breast_cancer_nb)
describe(breast_cancer_nb)

#summary statistics by class
describe(filter(breast_cancer_train, class_char == 'benign'))
describe(filter(breast_cancer_train, class_char == 'malignant'))

#visualise discretised distributions
breast_cancer_nb_v <- breast_cancer_nb

# assign id field for visualisations
breast_cancer_nb_v$id <- 1:nrow(breast_cancer_nb_v)

breast_cancer.m.d <- melt(select(breast_cancer_nb_v, -class_char),
                        id.var="id")

#visualise
ggplot(data = breast_cancer.m.d, 
       aes(x=as.numeric(value)),
       fill=class) + 
  geom_bar(aes(fill=factor(variable))) +
  facet_wrap("variable") +
  ggtitle("UCI Breast Cancer Dataset Distributions - DISCRETISED")


#split dataset into TRAIN and TEST

#remove variables not required
breast_cancer_train <- sample_frac(breast_cancer_nb, 0.8)
myrows <- as.numeric(rownames(breast_cancer_train)) 
breast_cancer_test <- breast_cancer_nb[-myrows,]

#inspect split datasets
str(breast_cancer_train)
str(breast_cancer_test)


#train naiveBayes model using e1071 package
nb_discretise <- naiveBayes(as.factor(class_char)~.,
                         data=breast_cancer_train,
                         laplace = 0)

#view model summary
nb_discretise
summary(nb_discretise)

#predictions
#prediction function
myPredictor <-function(nb_object, newdata) {
  predict(nb_object,
          newdata=newdata[,-10])  
}

myPredict_dc <- myPredictor(nb_discretise,
                         breast_cancer_test)


#confusion matrix
table(prediction=myPredict_dc,
      true=breast_cancer_test$class_char)

#accuracy
mean(myPredict_dc == breast_cancer_test$class_char)



#10 fold cross-validation
breast_cancer_dc <- tune(naiveBayes,
                           as.factor(class_char)~.,
                           data = breast_cancer_train,
                           predict.func = myPredictor
)

#review Cross Validation output
breast_cancer_dc$best.parameters
breast_cancer_dc$best.performance
breast_cancer_dc$performances
breast_cancer_dc$best.model


#use cross validated best model
nb_discrete_xv <- breast_cancer_dc$best.model

#run against test dataset
#prediction object
myPredict_dc_xv <- myPredictor(nb_discrete_xv,
                            breast_cancer_test)


#confusion matrix
table(prediction=myPredict_dc_xv,
      true=breast_cancer_test$class_char)

#accuracy
mean(myPredict_dc_xv == breast_cancer_test$class_char)
