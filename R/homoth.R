homoth <- function (form=rbind(c(0,1),c(1,1)),centre=c(0,1),scale=0.5)
{
# center on 0
form<-form-rbind(centre,centre)
form<-form*scale
form<-form+rbind(centre,centre)
form
}

