## Loading packages
library("data.table")
library("D3partitionR")

## Reading data
  titanic_data = fread("train.csv")

##Agregating data to have unique sequence for the 4 variables
var_names=c('Sex','Embarked','Pclass','Survived')
data_plot=titanic_data[,.N,by=var_names]
data_plot[,(var_names):=lapply(var_names,function(x){data_plot[[x]]=paste0(x,' ',data_plot[[x]])
})]

## Plotting the chart
library("magrittr")
plot <-D3partitionR() %>%
  add_data(data_plot,count = 'N',steps=c('Sex','Embarked','Pclass','Survived')) %>%
  add_title('Titanic') %>%
  plot()
