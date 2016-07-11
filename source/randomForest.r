args <- commandArgs(trailingOnly = TRUE)

training = args[1]
testing = args[2]
writing = args[3]
flag = as.numeric(args[4])

suppressMessages(require('randomForest'))
t = read.table(training,header = TRUE)
t1 = read.table(testing,header = TRUE)

if (flag == 0){
  x = cbind(t[1],t[3],t[4],t[5])
  x1 = cbind(t1[1],t1[3],t1[4],t1[5])
}
if (flag > 0){
  x = cbind(t[1],t[2],t[3],t[4],t[5])
  x1 = cbind(t1[1],t1[2],t1[3],t1[4],t1[5])
}

y = t$status
r = suppressWarnings(randomForest(x,y))
p = predict(r,x1)

out = p #cbind(p,t1$status)
write.table(out,file = writing,row.names = FALSE,col.names=FALSE,append=FALSE,sep='\t')
