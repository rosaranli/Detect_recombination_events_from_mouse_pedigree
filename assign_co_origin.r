##purpose: for each CO event, assign parental origin. 
source("/well/myers/ranli/data/MouseSequencing_bam/bin/Phase_F4.R")
co=read.table("/well/myers/ranli/data/MouseSequencing_bam/F45/VCFfiles/Intersect_MGP/hmm/co/co_unique_f5_overlap_with_origin.txt",as.is=T,header=T) #read CO file
origin=rep(NA,nrow(co))

files=list.files(pattern="geno_SNP10000_",path="/well/myers/ranli/data/MouseSequencing_bam/F45/VCFfiles/Intersect_MGP/hmm") ##sparsed version of HMM results 

for(j in 1:length(files)){
	file=files[j]
	print(j)
	geno=read.table(paste0("/well/myers/ranli/data/MouseSequencing_bam/F45/VCFfiles/Intersect_MGP/hmm/",file),as.is=T,header=T)
	names=colnames(geno)
	names1=substr(names,2,20)
	samples=unlist(strsplit(file,"_"))
	if(nrow(geno)==7){
                geno=geno[-3,]
		samples=samples[-5]
        }
	if(samples[5]=="168"){
		samples=samples[-5]
		geno=geno[-3,]
	}
	if(samples[length(samples)]=="chr6" & (sum(diff(as.numeric(geno[1,]))!=0)>6 || sum(diff(as.numeric(geno[2,]))!=0)>6)){
		names=colnames(geno)
		names1=substr(names,2,20)
		index=which(as.numeric(names1)<37000000 | as.numeric(names1)>56000000)
		geno=geno[,index]
		names1=names1[index]
	}
	chr=samples[length(samples)]
	samples1=samples[-c(1,2,length(samples))]
	
	pro=makeposterior(geno[,1:ncol(geno)],rep(0.0000000000000125,ncol(geno)-1),rep(0.5,ncol(geno)),errorprob=0.000000001)
	a=pro[[1]]
	f_index=which(a[1,]>0.4)
	m_index=which(a[2,]>0.4)
	share_index=f_index[f_index %in% m_index]
	f_index=f_index[!f_index %in% share_index]
	m_index=m_index[!m_index %in% share_index]
	weird=which(diff(f_index)==1)
	if(length(weird)==2){
				f_index1=f_index[-c(weird[1],weird[2]+1)]
				f_index=c(f_index1,f_index[weird[2]])		
	}
	weird=which(diff(m_index)==1)
	if(length(weird)==2){
		m_index1=m_index[-c(weird[1],weird[2]+1)]
		m_index=c(m_index1,m_index[weird[2]])
	}
	co_index=c()
	for(i in 3:length(samples1)){
		sample=samples1[i]
		co_index1=which(co[,1]==chr & co[,58]==sample)
		co_index=c(co_index,co_index1)
		co_index1=which(co[,1]==chr & co[,58]==paste0(sample,"_replacement"))
		co_index=c(co_index,co_index1)
		co_index1=which(co[,1]==chr & co[,58]==paste0(sample,"old"))
                co_index=c(co_index,co_index1)
	}
	if(length(f_index)+length(m_index)+length(share_index)!=length(co_index)){
		weird_file=c(weird_file,file)
		print(paste0("weird file ",file," total ",length(weird_file)))
	#	next
	}
	co_sub=co[co_index,]
	f_index=unique(f_index)
	m_index=unique(m_index)
	if(length(f_index)){
		for(i in 1:length(f_index)){
			co1=names1[f_index[i]-1]
			co2=names1[f_index[i]]
			co_sub_index=which(co_sub[,2]>=co1 & co_sub[,3]<=co2)
			origin[co_index[co_sub_index]]="f"
		}
	}
	if(length(m_index)){
		for(i in 1:length(m_index)){
			co1=names1[m_index[i]-1]
			co2=names1[m_index[i]]
			co_sub_index=which(co_sub[,2]>=co1 & co_sub[,3]<=co2)
			origin[co_index[co_sub_index]]="m"
		}
	}
	if(length(share_index)){
		for(i in 1:length(share_index)){
			co1=names1[share_index[i]-1]
			co2=names1[share_index[i]]
			co_sub_index=which(co_sub[,2]>=co1 & co_sub[,3]<=co2)
			origin[co_index[co_sub_index]]="unknown"
		}
	}
}

info4=cbind(co,origin)

write.table(info4,file="/well/myers/ranli/data/MouseSequencing_bam/F45/VCFfiles/Intersect_MGP/hmm/co/co_unique_f5_overlap_with_origin_simple_and_hmm.txt",row.names=F,quote=F)
