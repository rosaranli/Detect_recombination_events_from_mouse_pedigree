#author:Ran Li  rosa.ranli@gmail.com
#Purpose:hmm algorithm to smooth the genetic backgroud of F2, F4, F5 mice
#Usage:Rscript find_hmm_pass_allow_small_chunk.r probability_for_each_geno_type_from_GATK_call.file output small_trunk_threshold
args<-commandArgs(trailingOnly=TRUE)
#Get the HMM path accordington to posterior probability
find_route<-function(prob_matrix){
	state_seq<-unlist(apply(prob_matrix,1,function(v){return(which.max(v)-1)}))
	return(state_seq)
}

##########################################################
get_info<-function(state_seq){
	index<-1 #position
	count<-1 #number of successive states
	state<-state_seq[1]
	for(i in 2:length(state_seq)){
		if(state_seq[i]==state_seq[i-1]){
			cur_pos=length(index)
			count[cur_pos]=count[cur_pos]+1
		}
		else{
			index=c(index,i)
			count=c(count,1)
			state<-c(state,state_seq[i])
		
		}
	}
	return(cbind(index,state,count))
}

get_bg<-function(info,thres){
	count.less<-which(info[,3]<thres)#find the trunk less than thres
	num.info<-nrow(info)
	num.less<-length(count.less)
	if(!num.less){ #test if the length of the count.less is 0
		return(rep(info[,2],info[,3]))
	}
	if(count.less[1]==1){#if the first trunk is less than threshold
		if(num.less==1){#if only the first trunk is less than threshold then return
			return(rep(info[,2],info[,3]))
		}
		else{#else remove the first trunck
			count.less<-count.less[-1]
			num.less<-num.less-1
		}
	}
	
	if(num.less==1){# if the number of small trunk is 1
		if(count.less[1]==num.info){# if the small trunk is the last trunk then return
			return(rep(info[,2],info[,3]))
		}
		###add: if the states of previous and the next trunk are different, then return
		if(info[count.less[1]-1,2]!=info[count.less[1]+1,2]){
			return(rep(info[,2],info[,3]))
		}
		info[count.less[1]-1,3]<-info[count.less[1]-1,3]+info[count.less[1],3]
		info<-info[-count.less[1],]
		return(rep(info[,2],info[,3]))
		#return(info)
	}
	# if the number of the small trunk is greater than 2
	aux.stack<-info[1:(count.less[1]-1),,drop=F]#this auxiliary stack is used to record the process of  merging small trunk
	stack.size<-nrow(aux.stack)
	for(i in 1:(num.less-1)){
		if(info[count.less[i]-1,2]!=info[count.less[i]+1,2]){#if the states of previous trunk and next trunk of count.less[i]th trunk are different,then don't merge and push the current trunk into stack
			aux.stack<-rbind(aux.stack,info[count.less[i],])
			stack.size<-nrow(aux.stack)#update stack.size
		}
		else{#else
			aux.stack[stack.size,3]<-aux.stack[stack.size,3]+info[count.less[i],3]#update the counts in current state
		}		
		if(count.less[i+1]!=count.less[i]+1){#if the next small trunk is not the next trunk
			aux.stack<-rbind(aux.stack,info[(count.less[i]+1):(count.less[i+1]-1),])#update stack
			stack.size<-nrow(aux.stack)#update stack.size
		}
	}
	if(count.less[num.less]==num.info){#if the last small trunk is the last trunk then push the last trunk into stack and return
		aux.stack<-rbind(aux.stack,info[num.info,])
		return(rep(aux.stack[,2],aux.stack[,3]))
	}
	#else if not
	if(info[count.less[num.less]-1,2]!=info[count.less[num.less]+1,2]){#if the states of previous and next trunk of current trunk are different, then not merge and push the rest trunks into the stack and return
		aux.stack<-rbind(aux.stack,info[count.less[num.less]:num.info,])
		return(rep(aux.stack[,2],aux.stack[,3]))
	}
	aux.stack[stack.size,3]<-aux.stack[stack.size,3]+info[count.less[num.less],3]#merge the final small trunk
	aux.stack<-rbind(aux.stack,info[(count.less[num.less]+1):num.info,])
	return(rep(aux.stack[,2],aux.stack[,3]))
	#return(aux.stack)
}

input=read.table(args[1])
output=args[2]
thres=as.numeric(args[3])
state_seq<-find_route(input[,3:5])
info<-get_info(state_seq)
bg<-get_bg(info,thres)
write.table(cbind(input[,1:2],bg),file=output,row.names=F,col.names=F,quote=F)

