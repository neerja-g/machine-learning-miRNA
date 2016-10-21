##calculating sensitivity and specificity 
## (c) Neerja Garikipati 2016

GSEtable=(GSE61741table) #read in table

GSEmatrix=as.matrix(GSEtable) #assemble matrix


N=nrow(GSEmatrix)
M=ncol(GSEmatrix)
num_disease=73
num_normal= M-num_disease
sens=rep(0,N) #row vector of N 0s
spec=rep(0,N) #row vector of N 0s
pos_thresh=rep(0,N) #determined seperately for each miRNA

maxG=numeric(N) 
GSEmatrix_norm=matrix(data=NA, nrow=N, ncol=M, byrow=TRUE)

for(i in 1:N){
  maxG[i]=max(GSEmatrix[i,]) #maximum expression value of each miRNA (rows)
}

for(i in 1:N){
  for(j in 1:M){
    GSEmatrix_norm[i,j]=GSEmatrix[i,j]/maxG[i] #normalising between 0 and 1 for feature vector
  }
}


M_train_dis = 35 #Training set of disease samples
M_test_dis = num_disease - M_train_dis #Test set of disease samples
M_train_norm = 45 #Training set of normal samples
M_test_norm = num_normal - M_train_norm #Test set of normal samples
M_train = M_train_dis + M_train_norm
M_test = M_test_dis + M_test_norm

GSEmatrix_train=matrix(data=NA, nrow=N, ncol=M_train, byrow=TRUE) #define the training matrix
GSEmatrix_test=matrix(data=NA, nrow=N, ncol=M_test, byrow=TRUE) #define the test matrix

#Copy data into training set
for(i in 1:N){
  for(j in 1:M_train_dis){
    GSEmatrix_train[i,j] = GSEmatrix_norm[i,j]
  }
  for(j in M_train_dis+1:M_train_norm){
    j1 = num_disease + j - M_train_dis
    GSEmatrix_train[i,j] = GSEmatrix_norm[i,j1]
  }
}

#Copy data into testing set
for(i in 1:N){
  for(j in 1:M_test_dis){
    j1 = M_train_dis + j
    GSEmatrix_test[i,j] = GSEmatrix_norm[i,j1]
  }
  for(j in M_test_dis+1:M_test_norm){
    j1 = num_disease + M_train_norm + j - M_test_dis
    GSEmatrix_test[i,j] = GSEmatrix_norm[i,j1]
  }
}
#determining pos_thresh individually for each miRNA as the value maximising sens+spec
for(i in 1:N){
  pos_thresh_test = 0
  prev_max = 0.0
  while (pos_thresh_test < max(GSEmatrix_train[i,])){
    for(j in 1:M_train_dis){
      if(GSEmatrix_train[i,j] > pos_thresh_test){
        sens[i]= sens[i]+1
      }
    }
    for(j in M_train_dis+1:M_train_norm){
      if(GSEmatrix_train[i,j] <= pos_thresh_test){
        spec[i]=spec[i]+1
      }
    }
    sens[i]=sens[i]/M_train_dis
    spec[i]=spec[i]/M_train_norm
    if (sens[i] + spec[i] > prev_max){
      prev_max = sens[i] + spec[i]
      pos_thresh[i] = pos_thresh_test
    }
    pos_thresh_test = pos_thresh_test + 0.01
  }
#  print(paste('i=',i,"pos_thresh=",pos_thresh[i]))
}
#recompute sens and spec for each miRNA using pos_thresh for that miRNA
sens = rep(0,N) #vector of N 0s
spec = rep(0,N) #vector of N 0s
for(i in 1:N){
  for(j in 1:M_train_dis){
    if(GSEmatrix_train[i,j] > pos_thresh[i]){
      sens[i]= sens[i]+1   #sens = disease
    }
  }
  for(j in M_train_dis+1:M_train_norm){
    if(GSEmatrix_train[i,j] <= pos_thresh[i]){
      spec[i]=spec[i]+1  #spec = normal 
    }
  }
  sens[i]=sens[i]/M_train_dis
  spec[i]=spec[i]/M_train_norm
  print(paste('i=',i,"sens=",sens[i],"spec=",spec[i]))
}
index_sens = order(-sens) #descending order of sens
index_spec = order(-spec) #descending order of spec
for(i in 1:N){
  print(paste("i=",index_sens[i],"sens=",sens[index_sens[i]],"i=",index_spec[i],"spec=",spec[index_spec[i]]))
}


#Check disease samples in test set for 10 most sensitive miRNAs
for(i in 1:10){
  test_dis_sens = 0
  test_dis_spec = 0
  for(j in 1:M_test_dis){
    if(GSEmatrix_test[index_sens[i],j] > pos_thresh[index_sens[i]]){
      test_dis_sens = test_dis_sens + 1 #miRNA sensitive on disease sample
    }      
    if(GSEmatrix_test[index_spec[i],j] <= pos_thresh[index_spec[i]]){
      test_dis_spec = test_dis_spec + 1 #miRNA specific on disease sample
    }
    #    print(paste("By 10 most sens miRNA: Sample #, miRNA #",j, index_sens[i], GSEmatrix_test[index_sens[i],j]))
    #    print(paste("By 10 most spec miRNA: Sample #, miRNA #",j, index_spec[i], GSEmatrix_test[index_spec[i],j]))
  }
  test_dis_sens = test_dis_sens/M_test_dis
  test_dis_spec = test_dis_spec/M_test_dis
  print(paste("test_dis_sens, miRNA#:", test_dis_sens, index_sens[i], "test_dis_spec, miRNA#:", test_dis_spec, index_spec[i]))
}


#Check normal samples in test set for 10 most sens+spec miRNAs
for(i in 1:10){
  test_norm_sens = 0
  test_norm_spec = 0
  for(j in M_test_dis+1:M_test_norm){
    if(GSEmatrix_test[index_sens[i],j] > pos_thresh[index_sens[i]]){
      test_norm_sens = test_norm_sens + 1 #miRNA sensitive on normal sample
    }
    if(GSEmatrix_test[index_spec[i],j] <= pos_thresh[index_spec[i]]){
      test_norm_spec = test_norm_spec + 1 #miRNA specific on normal sample
    }
    
    
    #    print(paste("By 10 most sens miRNA: Sample #, miRNA #",j, index_sens[i], GSEmatrix_test[index_sens[i],j]))
    #    print(paste("By 10 most spec miRNA: Sample #, miRNA #",j, index_spec[i], GSEmatrix_test[index_spec[i],j]))
  }
  test_norm_sens = test_norm_sens/M_test_norm
  test_norm_spec = test_norm_spec/M_test_norm
  print(paste("test_norm_sens, miRNA#", test_norm_sens, index_sens[i], "test_norm_spec, miRNA#", test_norm_spec, index_spec[i]))
}


    
  
