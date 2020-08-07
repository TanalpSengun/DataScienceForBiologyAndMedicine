#Author: Tanalp Sengun
#Course: Data Science for Biology and Medicine
#Name of the Homework: Revealing Liver Cancer Factors from the Data

setwd("~/Desktop/clinical.project-TCGA-LIHC.2019-03-16-2")
# make sure that the finding is repeatable
set.seed(7)

#Libraries I have to upload
# load the library
library("mlbench")
library("caret")

#library(data.table)

library(glmnet)
#First I have to read the data from the txt
my_clinical_data <- read.table(file = "clinical.tsv", sep = "\t", header=TRUE,
                               quote = "\"", na.strings = "NA", fill = T, stringsAsFactors = T)
my_clinical_data2 <- read.table(file = "clinical.tsv", sep = "\t", header=TRUE,
                                quote = "\"", na.strings = "NA", fill = T, stringsAsFactors = T)

#After I look at my data I can see that therapeutic_agents,
#treatment_intent_type,treatment_or_therapy,days_to_recurrence
#days_to_last_known_disease_status are empty so they dont have importance
na_list=list(17,18,26,27,28)

#After I look at my data I can see that site_of_resection_or_biopsy,
#prior_malignancy,progression_or_recurrence,tumor_grade,morphology
#last_known_disease_status,classification_of_tumor,project_id 
#are all same so dont have a significance. Also ID's are not important
insignificance_list=list(1,2,3,9,10,15,19,20,22,23,24)

#columns to be deleted
columns_to_be_deleted_list=list(17,18,26,27,28,1,2,3,9,10,15,19,20,22,23,24)
columns_to_be_deleted_list<-as.numeric(columns_to_be_deleted_list)

#I substract the unnecessary data from my data so that it will increase the efficiency and the speed
my_clinical_data<-my_clinical_data[, -columns_to_be_deleted_list]
my_clinical_data2<-my_clinical_data2[, -columns_to_be_deleted_list]
#Now I have less feature that I can work on.

#I will turn every string into numeric values so that my algorthyms can work

my_clinical_data[,1]<-ifelse(my_clinical_data[,1]=='male', 1,0)

#afterwards males are 1, females are 0.
hist(my_clinical_data[,1], 
     main="Gender of Patients", 
     xlab="Gender", 
     border="blue", 
     col="green",
     xlim=c(0,1),
     las=2, 
     breaks=4)

library(plyr)

number_age = as.factor(revalue(my_clinical_data[,2], c("--"=1990)))
age_histogram = as.numeric(as.matrix(number_age))
#afterwards males are 1, females are 0.
hist(age_histogram, 
     main="Birth Year of Patients", 
     xlab="Birth_year", 
     border="blue", 
     col="green",
     xlim=c(1918,2000),
     las=1, 
     breaks=15)

#binarys

# Ggplot2 library
library(ggplot2)
age_at_diagnosis = as.factor(revalue(my_clinical_data[,8], c("--"=0)))
age_at_diagnosis = as.numeric(as.matrix(age_at_diagnosis))
hist(age_at_diagnosis, 
     main="age_at_diagnosis", 
     xlab="Birth_year", 
     border="blue", 
     col="green",
     xlim=c(10000,30000),
     las=1, 
     breaks=100)



#afterwards alives are 1, deads are 0.
my_clinical_data[,9]<-ifelse(my_clinical_data[,9]=='alive', 1,0)
vital_status<-my_clinical_data[,9]
Dead_Alive <- c(length(which(vital_status==1)),length(which(vital_status==0)))
lbls <- c("Alive", "Dead")
pie(Dead_Alive,labels = lbls, main="Dead or Alive")

#multiples

my_clinical_data[,4]<-ifelse(my_clinical_data[,4]=='not', 1, my_clinical_data[,4])
#afterwards hispanic or latinos are 1, not hispanic or latinos are 2, not reporteds are 3.

my_clinical_data[,3]<-ifelse(my_clinical_data[,3]=='indian', 1, my_clinical_data[,3])
#whites: 5       asian: 2     blacks: 3   indians:  1  not reporteds: 4

my_clinical_data[,7]<-ifelse(my_clinical_data[,7]=='iiib', 1, my_clinical_data[,7])
#stage 1: 2   #stage 2: 3   #stage 3: 4       #stage 3a:5    #stage 3b: 6   #stage 3c:7    #stage 4:8     #stage 4a:9    #stage 4b:10   #not reported:1        

my_clinical_data[,6]<-ifelse(my_clinical_data[,6]=='cholangio', 1, my_clinical_data[,6])

#after this I have completed to turn every string into number.


len=length(my_clinical_data[,4])
#more than 2 strings

#for(i in 1:len) {
#     if(my_clinical_data[i,4]=='not hispanic or latino'){
#      my_clinical_data[i,4]=1 }
#   else if(my_clinical_data[i,4]=='hispanic or latino'){
#    my_clinical_data[i,4]=2}
#else  {my_clinical_data[i,4]<=3}   
#}
my_clinical_data$days_to_death<- revalue(my_clinical_data$days_to_death, c("--"=0))
my_clinical_data$days_to_last_follow_up<- revalue(my_clinical_data$days_to_last_follow_up, c("--"=0))
my_clinical_data$year_of_death<- revalue(my_clinical_data$year_of_death, c("--"=0))
my_clinical_data$age_at_diagnosis<- revalue(my_clinical_data$age_at_diagnosis, c("--"=-1))
my_clinical_data$year_of_birth<-revalue(my_clinical_data$year_of_birth, c("--"=-1))

#to check the dependencies of variables and the survival I split them

all_other_factors<-my_clinical_data[,-9]

#this is the infromation in the end
str(my_clinical_data)


all_other_factors_matrix<-matrix(data=as.numeric(unlist(all_other_factors)),nrow = 377,ncol = 11 )
#this is the correlations of the factors with the vital status
res <- cor(all_other_factors_matrix,all_other_factors_matrix,method = c("pearson"))
round(res, 2)

res2 <- cor(all_other_factors_matrix,my_clinical_data$vital_status,method = c("pearson"))
round(res2, 2)

install.packages("Hmisc")
library("Hmisc")
res2 <- rcorr(as.matrix(all_other_factors_matrix))
res2
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P

library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#this is the correlation plot between the features.


all_other_factors_matrix<-matrix(data=as.numeric(unlist(all_other_factors)),nrow = 377,ncol = 11 )

Y<- as.matrix(vital_status)
lasso.fit<- glmnet(x=all_other_factors_matrix, y=Y, family = "binomial", alpha = 1)
plot(lasso.fit, xvar = "lambda",label = TRUE)

#In statistics and machine learning, lasso (least absolute shrinkage and selection operator; 
#also Lasso or LASSO) is a regression analysis method that performs both variable selection and 
#regularization in order 
#to enhance the prediction accuracy and interpretability of the statistical model it produces.

#lasso analizi masa??st??nde

#rm(my_clinical_data)

fisher_1<-fisher.test(all_other_factors_matrix[,1],vital_status)
fisher_2<-fisher.test(all_other_factors_matrix[,2],vital_status,simulate.p.value=TRUE,B=1e5)
fisher_3<-fisher.test(all_other_factors_matrix[,3],vital_status)
fisher_4<-fisher.test(all_other_factors_matrix[,4],vital_status)
fisher_5<-fisher.test(all_other_factors_matrix[,5],vital_status,simulate.p.value=TRUE,B=1e5)
fisher_6<-fisher.test(all_other_factors_matrix[,6],vital_status)
fisher_7<-fisher.test(all_other_factors_matrix[,7],vital_status,simulate.p.value=TRUE,B=1e5)
fisher_8<-fisher.test(all_other_factors_matrix[,8],vital_status,simulate.p.value=TRUE,B=1e5)
fisher_9<-fisher.test(all_other_factors_matrix[,9],vital_status)
fisher_10<-fisher.test(all_other_factors_matrix[,10],vital_status,simulate.p.value=TRUE,B=1e5)
fisher_11<-fisher.test(all_other_factors_matrix[,11],vital_status)


list_feature=list(fisher_1[1]<0.05,fisher_2[1]<0.05,fisher_3[1]<0.05,fisher_4[1]<0.05,fisher_5[1]<0.05,fisher_6[1]<0.05,fisher_7[1]<0.05
                  ,fisher_8[1]<0.05,fisher_9[1]<0.05,fisher_10[1]<0.05,fisher_11[1]<0.05)

list_feature

#2,3,5,7,9,11
#age,race,tumor_stage

#from my study I found that the most important 3 features for increasing the survival rate are age,race,
#tumor stage.

#Now I will try something different, I will not only think that survival is important I will also count the 
#days so that maybe living longer can be explained differntly with the data I have if you have cancer.


the_day_until_death <- my_clinical_data[,10]
the_day_last_follow_up <- my_clinical_data[,12]

#Create function to identify all columns that need repair

total_day_alive <- as.numeric(unlist(the_day_last_follow_up))
total_day_dead <- as.numeric(unlist(the_day_until_death))

total_all <- total_day_alive + total_day_dead
total_all<-cbind(vital_status,total_all)
the_factors_for_the_second <- my_clinical_data[,c(1,2,3,4,6,7,8)]

correlation_vector = vector()
num_the_second_factors = matrix(data=as.numeric(unlist(the_factors_for_the_second)),nrow = 377,ncol = 7 )

for (i in 1:7){
  correlation_vector[i] <- cor(num_the_second_factors[,i],total_all)
}

correlation_vector

f_1<-fisher.test(num_the_second_factors[,1],vital_status,simulate.p.value=TRUE,B=1e5)
f_2<-fisher.test(num_the_second_factors[,2],vital_status,simulate.p.value=TRUE,B=1e3)
f_3<-fisher.test(num_the_second_factors[,3],vital_status,simulate.p.value=TRUE,B=1e5)
f_4<-fisher.test(num_the_second_factors[,4],vital_status,simulate.p.value=TRUE,B=1e5)
f_5<-fisher.test(num_the_second_factors[,5],vital_status,simulate.p.value=TRUE,B=1e5)
f_6<-fisher.test(num_the_second_factors[,6],vital_status,simulate.p.value=TRUE,B=1e5)
f_7<-fisher.test(num_the_second_factors[,7],vital_status,simulate.p.value=TRUE,B=1e3)
fishers<-cbind(f_1[1],f_2[1],f_3[1],f_4[1],f_5[1],f_6[1],f_7[1])
fishers
#now I will print what are the labels smaller than p<0.05 which are significant
indices=vector()
num=0
for (i in 1:7){
  if(fishers[,i]<0.05){
    num= num+1
    indices[num]<-i
  }
  
}
#this will give the labels of the datas they are significantly important.
names(the_factors_for_the_second[indices])

#I will use one another method finally 
names(my_clinical_data2)   
my_clinical_data2
gender          <-as.factor(my_clinical_data2[,1])
year_of_birth   <-as.factor(my_clinical_data2[,2])                    
race            <-as.factor(my_clinical_data2[,3]) 
ethnicity             <-as.factor(my_clinical_data2[,4]) 
year_of_death         <-as.factor(my_clinical_data2[,5]) 
primary_diagnosis   <-as.factor(my_clinical_data2[,6]) 
tumor_stage<-as.factor(my_clinical_data2[,7]) 
age_at_diagnosis     <-as.factor(my_clinical_data2[,8]) 
days_to_death<-as.factor(my_clinical_data2[,10]) 
days_to_birth<-as.factor(my_clinical_data2[,11]) 
days_to_last_follow_up<-as.factor(my_clinical_data2[,12]) 

vital_status<-as.factor(my_clinical_data2[,9]) 


xfactors <- model.matrix(vital_status ~ gender + tumor_stage + primary_diagnosis + race)[, -1]
x        <- as.matrix(data.frame(xfactors))

# Note alpha=1 for lasso only and can blend with ridge penalty down to
# alpha=0 ridge only.
glmmod <- glmnet(x, y=as.factor(vital_status), alpha=1, family="binomial")

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")

glmnet(x = x, y = as.factor(vital_status), family = "binomial", alpha = 1) 

coef(glmmod)[,10]

vital_status_plot = matrix(data=as.numeric(unlist(vital_status)),nrow = 377,ncol = 1 )
cv.glmmod <- cv.glmnet(x, y=vital_status_plot, alpha=1)
plot(cv.glmmod)

#DO THE PCA !!!!!

all_other_factors_matrix

# log transform 
log.factors <- all_other_factors_matrix
they_dead <- as.matrix(vital_status)

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
cancer.pca <- prcomp(log.factors,
                 center = TRUE,
                 scale. = TRUE) 

print(cancer.pca)

plot(cancer.pca, type = "l")

summary(cancer.pca)

#finally my aim is to look at the mean of the different groups:
#died - long   #died - short     #alive - long  #alive - short   

final_labeling = my_clinical_data2[,cbind(9,10,12)]
final_labeling$days_to_death<- revalue(final_labeling$days_to_death, c("--"=NA))
final_labeling$days_to_last_follow_up<- revalue(final_labeling$days_to_last_follow_up, c("--"=NA))
mean_of_deceased = sum(as.numeric(as.matrix(na.omit(final_labeling$days_to_death))))/length(na.omit(final_labeling$days_to_death))
mean_of_alive = sum(as.numeric(as.matrix(na.omit(final_labeling$days_to_last_follow_up))))/length(na.omit(final_labeling$days_to_last_follow_up))

#Finally now I know how is the average
# 672.1364 deceased      #785.4807  alive mean

