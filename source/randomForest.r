suppressMessages(require('randomForest'))

args <- commandArgs(trailingOnly = TRUE)

trainingFile = paste(args[1],"data/trainingSet/",args[2], sep='')
testingFile = args[3]
writingFile = args[4]

t = read.table(trainingFile,header = TRUE)
t1 = read.table(testingFile,header = TRUE)

if (args[2] == 'genericAll.txt'){
  x = cbind(t[1],t[3],t[4],t[5])
  x1 = cbind(t1[1],t1[3],t1[4],t1[5])
}else{
  x = cbind(t[1],t[2],t[3],t[4],t[5])
  x1 = cbind(t1[1],t1[2],t1[3],t1[4],t1[5])
}

y = t$status
r = suppressWarnings(randomForest(x,y))
p = predict(r,x1)

out = p #cbind(p,t1$status)
write.table(out,file = writingFile,row.names = FALSE,col.names=FALSE,append=FALSE,sep='\t')
