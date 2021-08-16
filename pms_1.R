unlink(".RData")


options(scipen=999)

library(MASS);
.libPaths("/ibrix/home8/yongli/RCODE")
library(Matrix)
.libPaths("/ibrix/home8/yongli/RCODE"); 
library(lars);
.libPaths("/ibrix/home8/yongli/RCODE");
library(glmnet);
.libPaths("/ibrix/home8/yongli/RCODE");
library(mht);
.libPaths("/ibrix/home8/yongli/RCODE");
library(ncvreg)
.libPaths("/ibrix/home8/yongli/RCODE");
library(c060)

n.c1<-1;   

pen<-"lasso"

n.c<-12;  rou.1<-(0.6);   tiny<-0.30;
n.y<-250; n.p<-2000; devi<-1; 
nla<-100;      

sss_seq<-c(10,1000)  


columnmatrix<-matrix(rep(1:n.p,each=n.p),n.p,n.p);  rowmatrix<-matrix(rep(1:n.p,n.p),n.p,n.p)
Cov.Matrix<-matrix(0,n.p,n.p); Cov.Matrix<-(rou.1)^abs(columnmatrix-rowmatrix)
X <- mvrnorm(n.y, rep(0,n.p), Cov.Matrix);    er<-rnorm(n.y)


 
tauu<-sqrt(n.y/(2*log(n.p)))                 
                
trueind<-sample(1:n.p,n.c,replace=FALSE)


mincoeff<-2;
tinycoeff<-tiny*mincoeff

true.coeff<-c(rep(mincoeff,n.c-n.c1), rep(tinycoeff,n.c1))
print(true.coeff)




betaa<-rep(0,n.p); betaa[trueind]<-true.coeff;       
mu<-X%*%betaa;
signalnoise<-sqrt(sum(mu^2)/n.y)


Y<-mu+er; 


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################




pencv_f<-function(X,Y, pen)
        {
        
        pencv_fit<-cv.ncvreg(X, Y, family=c("gaussian"),   penalty=c(pen))
        pencv_model<-(ncvreg(X, Y, family=c("gaussian"), penalty=c(pen), lambda=pencv_fit$lambda.min)$beta)[-1]
        pencv_model
        }

################################################################################

pen_f<-function(X,Y, pen)
      {
      
      pen.beta.hat<-ncvreg(X, Y, family=c("gaussian"),penalty=pen)$beta[-1,]                 
      pen.candidate<-unique(as.matrix(t(pen.beta.hat!=0)));
      pen.candidate
      
      }
################################################################################



sel_ricc<-function(inc_m)
          {
          
          candidate<-inc_m==1 ; dm<-apply(candidate,1,sum);  rss<-rep(0,nrow(candidate) );  
          
          
          for (jj in 1:nrow(candidate)) 
                                      {
                                      if (dm[jj]==0) {rss[jj]<-deviance(lm(Y~1))} else{  rss[jj]<-deviance(lm (Y~X[,candidate[jj,]]))  }
                                      };  
          
          
          ricc<-n.y*log(rss/n.y)+2*(log(n.p)+2*log(log(n.p)))*dm;
          
          
          sel_ricc<- (1:n.p)[candidate[(ricc==min(ricc)),]]
          
          
          sel_ricc
          }

################################################################################
bolasso_f<-function(X,Y, sss)
          {
          
          la_seq<-ncvreg(X, Y, family=c("gaussian"),  nlambda=nla, penalty="lasso")$lambda 
          mod<-bolasso(X,Y, mu= la_seq , m=sss, probaseuil=1)
          bolasso_inc<-(mod$ind)[-1,]
          
          bolasso_inc_m <- unique(t(bolasso_inc))
          bolasso_inc_m
          
          }
          
          
################################################################################
pms_inc_f<-function(X,Y,pen,sss)
          {
          
          coeff<-(ncvreg(X, Y, family=c("gaussian"), penalty=pen, lambda=  sqrt(log(n.p)/n.y) )$beta)[-1]
           
          
          pms.candidate1<-NULL
          
          for (k2 in 1:(sss))
                              {                  
                              PM<-(matrix(rnorm(n.y*n.p),n.y,n.p))*tauu
                              Yhat<-(PM)%*%(coeff)
                              
                              XX<-NULL;YY<-NULL
                              XX<-(X+PM);
                              YY<-Y+Yhat                    
                              pen.beta.hat<-ncvreg(XX, YY, family=c("gaussian"),  nlambda=nla, penalty=pen)$beta[-1,]                 
                              pms.candidate2<-as.matrix(t(pen.beta.hat!=0));                                        
                              pms.candidate1<-rbind(pms.candidate1,pms.candidate2)
                                                  
                              PM<-NULL; XX<-NULL;YY<-NULL;m2<-NULL;
                              }
          pms_inc_m<-unique(pms.candidate1); 
          pms_inc_m
          
                                                            
          }
          
################################################################################

stabsel_f<- function(X,Y, sss)
            {
            
            
            xm<-stabpath(Y,X,steps=sss,weakness=1,family=c("gaussian"),nlambda=nla)
            cand<-t(xm[[2]]==1 )
            candidate<-unique(cand)
            candidate
            }
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################






pencv_sel_time<-(system.time(pencv_model<-pencv_f(X,Y,pen)    ))[1]
pencv_sel_time


pencv_sel<- sum((pencv_model!=0)==(betaa!=0))==n.p
pencv_p <- sum((pencv_model!=0)[betaa==0])
pencv_n <- n.c-sum((pencv_model!=0)[betaa!=0])
pencv_out<-c(pencv_sel,pencv_p, pencv_n,pencv_sel_time)
pencv_out

################################################################################
################################################################################



pen_inc_time<-(system.time(pen_inc_m<-pen_f(X,Y,pen)))[1]
pen_inc_time

pen_size<-nrow(pen_inc_m)
pen_size

pen_inc<-max(apply((t(pen_inc_m)==(betaa!=0)),2,sum))  ==n.p
pen_inc



pen_sel_time<-(system.time(pen_model<-sel_ricc(pen_inc_m)))[1]
pen_sel_time


pen_sel<-setequal(pen_model, trueind)
pen_p<-length(setdiff(pen_model, trueind))
pen_n<-length(setdiff(trueind,pen_model))




pen_out<-c(pen_size,pen_inc,pen_inc_time, pen_sel,pen_p, pen_n,pen_sel_time)
pen_out


################################################################################
################################################################################



bo_output<-NULL ; pms_output<-NULL; ss_output<-NULL;

for (sss in sss_seq)


                    {
                    #sss<-10
                    ss_inc_time<-(system.time(stabse_inc_m<-stabsel_f(X,Y,sss)))[1]
                    ss_inc_time
                    
                    ss_size<-nrow(stabse_inc_m)
                    ss_size
                    
                    ss_inc<-max(apply((t(stabse_inc_m)==(betaa!=0)),2,sum))  ==n.p
                    ss_inc
                    
                    ss_sel_time<-(system.time(stabse_model<-sel_ricc(stabse_inc_m)))[1]
                    ss_sel_time
                    
                    
                    ss_sel<-setequal(stabse_model, trueind)
                    ss_p<-length(setdiff(stabse_model, trueind))
                    ss_n<-length(setdiff(trueind,stabse_model))
                    
                    ss_out<-c(sss, ss_size,ss_inc,ss_inc_time, ss_sel,ss_p, ss_n,ss_sel_time)
                    ss_out
                    
                    ss_output<-c(ss_output,ss_out)
                    
                    ################################################################################
                    ################################################################################
                    
                                    
                    bo_inc_time<-(system.time(bolasso_inc_m<-bolasso_f(X,Y,sss)))[1]
                    bo_inc_time
                    
                    bo_size<-nrow(bolasso_inc_m)
                    bo_size
                    
                    bo_inc<-max(apply((t(bolasso_inc_m)==(betaa!=0)),2,sum))  ==n.p
                    bo_inc
                    
                    bo_sel_time<-(system.time(bolasso_model<-sel_ricc(bolasso_inc_m)))[1]
                    bo_sel_time
                    
                    
                    bo_sel<-setequal(bolasso_model, trueind)
                    bo_p<-length(setdiff(bolasso_model, trueind))
                    bo_n<-length(setdiff(trueind,bolasso_model))
                    
                    bo_out<-c(sss, bo_size,bo_inc,bo_inc_time, bo_sel,bo_p, bo_n,bo_sel_time)
                    bo_out
                    
                    bo_output<-c(bo_output,bo_out)

                    
                    ################################################################################
                    ################################################################################


                    pms_inc_time<-(system.time(pms_inc_m<-pms_inc_f(X,Y,pen,sss)))[1]
                    pms_inc_time
                                                      
                    pms_size<-nrow(pms_inc_m)
                    pms_size
                    pms_inc<-max(apply((t(pms_inc_m)==(betaa!=0)),2,sum))==n.p
                    pms_inc
                    
                    pms_sel_time<-(system.time(pms_model<-sel_ricc(pms_inc_m)))[1]
                    pms_sel_time
                    
                    
                    pms_sel<-setequal(pms_model, trueind)
                    pms_p<-length(setdiff(pms_model, trueind))
                    pms_n<-length(setdiff(trueind,pms_model))
                    
                    
                    
                    pms_out<-c(sss, pms_size,pms_inc,pms_inc_time, pms_sel,pms_p, pms_n,pms_sel_time)
                    pms_output<-c(pms_output,pms_out)

                    
                    }


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
mre<-c(n.y,n.p, n.c, n.c1, tauu, rou.1, mincoeff, "cv",pencv_out, "pen",pen_out, "ss",ss_output, "bo",bo_output,"pms",pms_output)
print(mre)



vvv<-round(abs(rnorm(1)*100000000))
filename<-paste(vvv,".",n.c1,"txt",sep="")
write.table(t(mre), file=filename, row.names=FALSE,col.names=FALSE)


                                      