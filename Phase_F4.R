####analyse a data matrix given information on the data (ourdata), and recombination probabilities (dists)
####prior on SNP frequencies is freq
####0,1,2 format and code missing data as "3"
####first two rows of "ourdata" give known genotype information on parents (missing parents all "3")
####dists should be one shorter than ourdata
####freqs must have same length as number of columns in ourdata

makeposterior=function(ourdata,dists,freq,errorprob=0.01,initforward=NA, initbackward=NA){
	nchildren=nrow(ourdata)-2
	assign("nchildren",nchildren,envir=.GlobalEnv)
	n=nchildren
	###now record a set of emission probabilities with same dimension?
	emitprob=array(dim=c(16,4^nchildren,4*(2+nchildren)))
	underlyingstate=emitprob[,,1:(2*nchildren)]  ####by ran. it has nothing to do with emitprob, just borrowed its format
	for(i in 0:1) for (j in 0:1) for(k in 0:1) for(l in 0:1){
		possmat=c(i,j)
		posspat=c(k,l)
		parentstate=8*i+4*j+2*k+l
		for(m in 1:4^nchildren){
			###transmission vector
			###a vector with 2n bits, first n from mother, next n from father
			u=getvec(nchildren,(m-1))	
			###which chromosome from mother for each child
			###first n bits received from mother
			maternal=u[1:nchildren]+1
			###which chromosome from father for each child
			###last n bits received from father
			paternal=u[(nchildren+1):(2*nchildren)]+1
			#use these vectors to index into possmat and posspat.
			#copy which allele is on resp mother and father chrs into the child chr
			#ourvals is therefore a bitstring of length 2n, at the moment a faithful representation of the txn vector, errors will be dealt with soon	
			ourvals=c(possmat[maternal],posspat[paternal])
			underlyingstate[(8*i+4*j+2*k+l+1),m,]=ourvals
		}  ####end of for m
	} #####end of for i,j,k,l
	#####underlying state element i,j,k reads binary (i-1) for genotype of two parents mother then father, binary (j-1) for inheritance vector
	#####(first maternal then paternal) gives out alleles (first maternal, then paternal)
	####emissions prob
	#temp prob gives the probability of emitting a 1 in the relevant position of the child bistring. 
	tempprob=underlyingstate  ####dim 16x4^(nchildren)x(2*nchildren)
	tempprob[underlyingstate==1]=1-errorprob
	tempprob[underlyingstate==0]=errorprob
	#convert that into an emission probability
	####I think plus 8 means the first 8 should represent parents
	emitprob[,,8+seq(1,4*nchildren,4)]=(1-tempprob[,,1:nchildren])*(1-tempprob[,,(nchildren+1):(2*nchildren)])  ###seq(9+4*(2+nchildren),4)
	emitprob[,,8+seq(2,4*nchildren,4)]=(1-tempprob[,,1:nchildren])*(tempprob[,,(nchildren+1):(2*nchildren)])+(tempprob[,,1:nchildren])*(1-tempprob[,,(nchildren+1):(2*nchildren)])
	emitprob[,,8+seq(3,4*nchildren,4)]=(tempprob[,,1:nchildren])*(tempprob[,,(nchildren+1):(2*nchildren)])
	emitprob[,,8+seq(4,4*nchildren,4)]=1
	probs=c(errorprob,1-errorprob)
	for(i in 0:1) for (j in 0:1) for(k in 0:1) for(l in 0:1){
 		curprobs=c(probs[i+1],probs[j+1],probs[k+1],probs[l+1])
curemits=c( (1-curprobs[1])*(1-curprobs[2]) ,
 (curprobs[1])*(1-curprobs[2])+(1-curprobs[1])*(curprobs[2])  ,
(curprobs[1])*(curprobs[2]) ,
1,
(1-curprobs[3])*(1-curprobs[4]) ,
 (curprobs[3])*(1-curprobs[4])+(1-curprobs[3])*(curprobs[4])  ,
(curprobs[3])*(curprobs[4]) ,
1
)

		possmat=c(i,j)
		posspat=c(k,l)

		parentstate=8*i+4*j+2*k+l+1
		for(m in 1:4^nchildren){
			emitprob[parentstate,m,1:8]=curemits      ######no matter what m is, parents part does not change
		}
	} ####end of i,j,k,l
	assign("emitprob",emitprob,envir=.GlobalEnv)
	######the probability of each child emitting each state, given an underlying state vector
	#####only allow single recombination events
	####just worry about children
	###no event matrix
	###multiply by prob of no event
	nstates=4^{nchildren}
	noevent=matrix(0,nrow=nstates,ncol=nstates)
	noevent[cbind(1:nstates,1:nstates)]=1  ###### assign the diagonal to 1
	####need 1 minus prob of no event times one over nchildren to multiply
	oneevent=noevent*0
	for(m in 1:4^nchildren){
      	for(n in 1:4^nchildren){
			u=getvec(nchildren,(m-1))
			v=getvec(nchildren,(n-1))
			ourdist=sum(abs(u-v))
			if(ourdist==1) oneevent[m,n]=1
		}
	}
	onematevent=noevent*0
	for(m in 1:4^nchildren){
      	for(n in 1:4^nchildren){
			u=getvec(nchildren,(m-1))
			v=getvec(nchildren,(n-1))
			ourdist=sum(abs(u-v)[1:nchildren])
			if(ourdist==1) onematevent[m,n]=1
		}
	}
	onepatevent=noevent*0
	for(m in 1:4^nchildren){
      	for(n in 1:4^nchildren){
			u=getvec(nchildren,(m-1))
			v=getvec(nchildren,(n-1))
			ourdist=sum(abs(u-v)[(nchildren+1):(2*nchildren)])
			if(ourdist==1) onepatevent[m,n]=1
		}
	}
	####not both
	temp=onematevent
	onematevent[oneevent!=1]=0
	onepatevent[oneevent!=1]=0
	####from m to n
	#######this is all preliminary work
	#####to step, calculate current child state terms, then multiply to get next child state terms
	#####also get prior parent state terms. mult by emissions probs
	#####for initial state probs have vector of constants, 1/4^nchildren.
	oneevent=onematevent+onepatevent
	assign("onematevent",onematevent,envir=.GlobalEnv)
	assign("onepatevent",onepatevent,envir=.GlobalEnv)
	assign("oneevent",oneevent,envir=.GlobalEnv)
	stateprob=array(dim=c(16,4^nchildren))
	nsites=ncol(ourdata)
	##posterior=array(dim=c(16,4^nchildren,nsites))
	forward=array(dim=c(16*3,4^nchildren,nsites))
	backward=array(dim=c(16,4^nchildren,nsites))
	constforward=vector(length=nsites)
	###look at recombination probability
	matprob=constforward
	patprob=matprob
	###forward
	#####assume have an n+2 by nsites matrix of genotypes
	#####also have an nsites vector of frequencies
	#####also have an nsites-1 vector of recombination probabilities
	####initialise stateprob for forward
	
	###parental state probability
	#freq is the allele frequency of the "1" allele 
	priorstateprob=matrix(nrow=16,ncol=nsites)
	for(i in 0:1) for(j in 0:1) for(k in 0:1) for(l in 0:1){
		total1=i+j+k+l
		total0=4-total1
		priorstateprob[8*i+4*j+2*k+l+1,]=freq^total1*(1-freq)^total0
	}
	assign("priorstateprob",priorstateprob,envir=.GlobalEnv)
	
	#######new distances
	dist=dists[1]
	newdata=ourdata[,1]
	pos=1
	if(is.na(initforward)) for(i in 1:16) for(j in 1:4^nchildren) stateprob[i,j]=1/4^nchildren/16
	else stateprob=initforward;
	pos=0
	const=0;
	const=1:nsites*0
	for(i in 1:nsites){
		pos=pos+1
		newdata=ourdata[,pos]
		#if(!pos %%100) print(pos)
		if(pos>1) dist=dists[(pos-1)]
		else dist=0.0;
		v=stepforward(stateprob,dist,newdata,pos)
		stateprob=v[1:16,]
		if(pos>1) const[pos]=const[pos-1]
		if(max(stateprob)<1e-5){
		const[pos]=const[pos]+log(max(stateprob))
		v=v/max(stateprob)
		stateprob=stateprob/max(stateprob)
		}
	forward[,,pos]=v
	}
	####need backward probs 
	####encouragingly fast, should do up to 4-5 children?

	dist=dists[length(dists)]
	newdata=ourdata[,ncol(ourdata)]
	pos=nsites
	for(i in 1:16) for(j in 1:4^nchildren) stateprob[i,j]=1
	backward[,,pos]=stateprob
	const2=1:nsites*0
	for(i in nsites:2){
		pos=pos-1
		if((pos+1)<=ncol(ourdata)) newdata=ourdata[,(pos+1)]
		else newdata=rep(4,nchildren+2)
		dist=dists[pos]
		v=stepbackward(stateprob,dist,newdata,pos)
		stateprob=v;
		const2[pos]=const2[pos+1]
		if(max(stateprob)<1e-5){
		const2[pos]=const2[pos]+log(max(stateprob))
		stateprob=stateprob/max(stateprob)
		}
		backward[,,pos]=stateprob
	}

	###make forward into posterior
	forward[1:16,,]=forward[1:16,,]*backward
	forward[17:32,,]=forward[17:32,,]*backward
	forward[33:48,,]=forward[33:48,,]*backward
	overallconst=const+const2
	overallsum=overallconst*0
	for(i in 1:nsites) overallsum[i]=sum(forward[1:16,,i])
	overallconst=overallconst+log(overallsum)
	#####same likelihood everywhere, a good thing
	for(i in 1:nsites) forward[,,i]=forward[,,i]/overallsum[i]
	###recombination probability
	matrecprob=overallconst*0
	patrecprob=matrecprob
	for(i in 1:nsites){
 		matrecprob[i]=sum(forward[17:32,,i])
		patrecprob[i]=sum(forward[33:48,,i])
	}
	###parentalposterior
	parentalposterior=matrix(nrow=16,ncol=nsites)
	for(i in 1:16) for(j in 1:nsites) parentalposterior[i,j]=sum(forward[i,,j])
	a=rbind(matrecprob,patrecprob,parentalposterior,overallconst)
	b=forward[1:16,,]
	return(list(a,b))
	###end of makeposterior
}

getvec=function(nchildren, childind){
	first=childind%%2
	childind=(childind-first)/2
	for(i in 2:(2*nchildren)){
		news=childind%%2
		childind=(childind-news)/2
		first=c(news,first)
	}
	return(first)
}

####newdata will use code 0 1 2 3 for genotypes and missing data
stepforward=function(stateprob,dist,newdata,pos,jump=T){
	nchildren=length(newdata)-2
	norec=exp(-(dist*nchildren*2))
	rec=(1-norec)/2/nchildren
	####probability of data given each state
	indices=seq(1,4*(nchildren+2),4)+newdata ##length nchildren+2
	curemit=emitprob[,,indices[1]]
	for(i in 2:length(indices)) curemit=curemit*emitprob[,,indices[i]]  ######times probability of genotype for each sample one by one
	####probability of going from each transmission vector state; transmission vector is independent of the state assignment of the parental states. Therefore colSums over the 16 parental states to one row with 4^nchildren columns  
	newprob=matrix(colSums(stateprob),nrow=1)
	####premultiply
	getprobnochange=norec*newprob
	getprobmatchange=rec*newprob%*%onematevent  ##(Ran)no matter what state it came from ,the likelihood of each state
	getprobpatchange=rec*newprob%*%onepatevent
	getprob=getprobnochange+getprobpatchange+getprobmatchange
	####fraction of likelihood from a change in state between i-1 and i
	matcont=getprobmatchange/getprob
	patcont=getprobpatchange/getprob
	newstateprob=curemit
	matprob=newstateprob
	patprob=matprob
	for(i in 1:16){
		#newstateprob on RHS has emission probabilities for this parental state and each  
		#for the parental states since each state is assumed to be independent of the previous state, the transition probability of going from one state to another is just the probability of the latter state in terms of  allele frequency	
		newstateprob[i,]=newstateprob[i,]*getprob*priorstateprob[i,pos]
		matprob[i,]=newstateprob[i,]*matcont
		patprob[i,]=newstateprob[i,]*patcont
	}
	stateprob=newstateprob
	if(jump==F) return(stateprob)
	else if(jump==T){
		v=array(dim=c(dim(stateprob)[1]*3,dim(stateprob)[2]))
		v[(1:dim(stateprob)[1]),]=stateprob
		v[((1:dim(stateprob)[1])+dim(stateprob)[1]),]=matprob
		v[((1:dim(stateprob)[1])+dim(stateprob)[1]*2),]=patprob
		return(v)
	}
} ### end of function stepforward

####newdata will use code 0 1 2 3 for genotypes and missing data
stepbackward=function(stateprob,dist,newdata,pos){
	nchildren=length(newdata)-2
	norec=exp(-(dist*nchildren*2))
	rec=(1-norec)/2/nchildren
	####probability of data given each state
	indices=seq(1,4*(nchildren+2),4)+newdata
	curemit=emitprob[,,indices[1]]
	for(i in 2:length(indices)) curemit=curemit*emitprob[,,indices[i]]
	####now premultiply by emission probability - for previous state
	stateprob=stateprob*curemit
	for(i in 1:16){
		stateprob[i,]=stateprob[i,]*priorstateprob[i,(pos+1)]
	}
	####probability of going to each state
	newprob=matrix(colSums(stateprob),nrow=1)
	####premultiply
	getprobnochange=norec*newprob
	getprobchange=rec*newprob%*%oneevent
	getprob=getprobnochange+getprobchange
	####fraction of likelihood from a change in state between i-1 and i
	##matcont=getprobmatchange/getprob
	##patcont=getprobpatchange/getprob
	newstateprob=stateprob*0
	for(i in 1:16){
		newstateprob[i,]=getprob
	}
	stateprob=newstateprob
	return(stateprob)
} #### the end of the function stepbackward 



