#author:rosa.ranli@gmail.com
#Purpose:to get the potential converted sites. Input are smoothed HMM results and genotypes

args<-commandArgs(trailingOnly=TRUE)
find_route<-function(prob_matrix){ 
	#state_seq<-apply(prob_matrix,1,function(v){return(which(v==0)[1]-1)})
	state_seq<-apply(prob_matrix,1,function(v){return(which(v==0)-1)})
	return(state_seq)
}

##########################################################

called=read.table(args[1],as.is=T,fill=T)
filter<-apply(called,1,function(m)if(any(is.na(m))) return (F) else return(T))
called<-called[filter,]
#filter2<-apply(raw_small[,3:5],1,function(m)if(any(m==0)) return (T) else return(F))
#raw_small<-raw_small[filter2,]
filter3<-apply(called[,3:5],1,function(m)if(length(which(m==0))==1) return (T) else return(F))
called<-called[filter3,]
called<-called[order(called[,2]),]
hmm_bg=read.table(args[2])
called_ems<-called[,3:5]
called_state_route<-find_route(called_ems)
full_info<-cbind(called[,1:2],called_state_route,hmm_bg[,2:3])
#full_info=cbind(called_state,hmm_bg[,3])
#full_info=cbind(called_state,hmm_bg[,2:3])
diff_index=which(full_info[,3]!=full_info[,5])
diff_info=full_info[diff_index,]
if(diff_info[1,1]!="chr6" || (diff_info[1,1]=="chr6" && sum(diff_info[,2]>37000000 & diff_info[,2]<56000000)/nrow(diff_info)<0.5)){
	write.table(diff_info,file=args[3],row.names=F,col.names=F,quote=F)
#write.table(full_info,file=args[4],row.names=F,col.names=F,quote=F)
}else{
	write.table(diff_info,file=paste0(args[3],"_original"),row.names=F,col.names=F,quote=F)
	diff_info2=diff_info[diff_info[,2]<37000000 | diff_info[,2]>56000000,]
	write.table(diff_info2,file=args[3],row.names=F,col.names=F,quote=F)
}
