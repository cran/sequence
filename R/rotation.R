rotation <- function (form,center,angle) 
{
npoints<-length(form[,1])
for (i in 1:npoints)
	{
	 form[i,]<-form[i,]-center
	}
matrot<-matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2)
rotated<-matrot %*% t(form) + cbind(center,center)
rotated<-data.matrix(t(rotated))
row.names(rotated)<-NULL
rotated
}

