###author: rosa.ranli@gmail.com

args<-commandArgs(trailingOnly=TRUE)
ori_file=read.table(args[1]) #this file shows the smoothed HMM results for each SNP
sample=(args[2])
out=args[3]
get_crossover<-function(a){
	pos<-as.vector(as.numeric(a[,3]))
	diff=diff(pos)
	co_start_index<-which(diff!=0)
	co_end_index=co_start_index+1
	#result=matrix(0,nrow=length(co_start_index)+1,ncol=5)
	result=matrix(0,nrow=length(co_start_index),ncol=7)
	#first column: chromosome number
	#column2: crossover(co) start position
	#column3: crossover(co) end position
	#column5: the number of snps from the upstreat crossover or the beginning of the chromosome
	#column4: the background upstream the crossover
	#result has length(co_start)_1 rows because although there is length(co_start) state change but there is plus 1 background state. I add one more row to record the information from the last crossover end point to the end of the chromosome, just as the first row
	chr=as.vector(a[1,1])
	result[,1]=chr
	result[,2]=a[co_start_index,2]
	result[,3]=a[co_end_index,2] 
	result[,4]=a[co_start_index,3]
	result[,5]=a[co_end_index,3]
	interval_up=c(0,co_start_index)
	interval_down=c(co_end_index,nrow(a))
	co_length_up=diff(interval_up)#get the number of snps for each bg
	co_length_down=diff(interval_down)
	result[,6]=co_length_up
	result[,7]=co_length_down
	
	return(result)
}
result=get_crossover(ori_file)
result=cbind(result,sample)
write.table(result,file=out,row.names=F,col.names=F,quote=F,append=T)




