read.old.seq <- function (file="NULL") 
{
seriseq<-list(NULL)
serie<-scan(file,what="character")
n<-length(serie)
message("\n Number of lines read: ",n,"\n")
listind<-c(NULL)
j<-0
seriseq[[1]][1]<-serie[1]
listind[1]<-1
for (i in 2:n)
	{
	 if(substring(serie[i],1,1)=="%") 
		{ 
		 j<-j+1
             # message("\nj=",j,"\n")
             seriseq<-c(seriseq,serie[i])
             listind[j+1]<-i
            }
      }
# message("\n listind \n\n")
print(listind)
nind<-length(listind)
print(nind)
for(i in 2:(nind+1))
    { 
	if(i==(nind+1)) m<-n-listind[i-1] else m<-listind[i]-1-listind[i-1]
      for (k in 1:m) seriseq[[i-1]][k+1]<-serie[listind[i-1]+k]
                        
    }

seriseq
}

