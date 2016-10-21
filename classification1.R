##classification algorithm for miRNA expression values
## (c) Neerja Garikipati 2016
GSEtable=(GSE61741table)

GSEmatrix=as.matrix(GSEtable)

N=nrow(GSEmatrix)
M=ncol(GSEmatrix)
num_disease=73
num_normal= M-num_disease

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
num_disease=73
num_normal= M-num_disease
M_train_dis = 35 #Training set of disease samples
M_test_dis = num_disease - M_train_dis #Test set of disease samples
M_train_norm = 45 #Training set of normal samples
M_test_norm = num_normal - M_train_norm #Test set of normal samples
M_train = M_train_dis + M_train_norm
M_test = M_test_dis + M_test_norm

GSEmatrix_train=matrix(data=NA, nrow=N, ncol=M_train, byrow=TRUE)
GSEmatrix_test=matrix(data=NA, nrow=N, ncol=M_test, byrow=TRUE)

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

theta=runif(N, -1.0, 1.0) #initizialise theta vector to randoms between 0 and 1

y_dis=rep(1, M_train_dis) #Training data for disease samples = 1
y_norm=rep(0, M_train_norm) #Training data for normal samples = 0

y=c(y_dis, y_norm) #vector of sample results

alpha=0.01 #learning rate. Good value =0.01
lambda=0.00 #penalization parameter to prevent overfitting to training data. Good value=0

theta_T=t(theta)

costFunc = 0.5*lambda*sum(theta*t(theta)) #initialise cost function to penalisation 

h=numeric(M_train)
for(j in 1:M_train){
  h[j] = 1/(1 + exp(-sum(theta*GSEmatrix_train[,j])))
  if(h[j] != 0 & h[j] != 1){
    costFunc = costFunc - (y[j]*log(h[j]) + (1 - y[j])*log(1 - h[j]))
  }
#  print(paste("j=", j, "h=", h[j], "cost function=", costFunc))
}

#Good value: costFunc <= 0.1
while(costFunc > 0.05){
  for(i in 1:N){
    costFunc_theta = lambda*theta[i]
    for(j in 1:M_train){
      costFunc_theta = costFunc_theta + (h[j] - y[j])*GSEmatrix_train[i,j]
    }
    theta[i] = theta[i] - alpha*costFunc_theta
  }
  costFunc = 0.5*lambda*sum(theta*t(theta))
  for(j in 1:M_train){
    h[j]=1/(1 + exp(-sum(theta*GSEmatrix_train[,j])))
    if(h[j] != 0 & h[j] != 1){
      costFunc = costFunc - (y[j]*log(h[j])+ (1 - y[j])*log(1 - h[j]))
    }
#    print(paste("j=", j, "h=", h[j], "cost function=", costFunc))
  }
  print(paste("cost function=", costFunc))
}

#x_trial = runif(N, 0.0, 1.0) # generate a random vector of N miRNA values from 0 to 1
#h_trial = 1/(1+exp(-sum(theta*t(x_trial))))

#Testing with samples from GSEmatrix_test
for(j in 1:M_test){
  x_test = GSEmatrix_test[,j]
  h_test = 1/(1+exp(-sum(theta*x_test)))
  print(paste("j=", j, "probability of positive outcome=", h_test))
}



