# Frozen-model evaluation; all fitted prediction preprocessing comes from set-a.
suppressPackageStartupMessages({library(tidymodels);library(data.table);library(ggplot2);library(knitr)})
source('checkpoint3_features.R')
data.table::setDTthreads(2)
set.seed(1707)
started <- Sys.time()
dir.create('output',showWarnings=FALSE)
save_table <- function(x,name) readr::write_csv(x,file.path('output',paste0('cp3_',name,'.csv')))
save_plot <- function(x,name,w=9,h=6) ggsave(file.path('output',paste0('cp3_',name,'.png')),x,width=w,height=h,dpi=170,bg='white')
fmt <- function(x,d=3) formatC(x,format='f',digits=d)
auc <- function(y,p) {
  n1<-sum(y==1); n0<-sum(y==0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)
}
boot_auc <- function(y,p,B=2000) {
  i0<-which(y==0);i1<-which(y==1)
  replicate(B,{i<-c(sample(i0,length(i0),TRUE),sample(i1,length(i1),TRUE));auc(y[i],p[i])})
}
cal <- function(y,p) {
  lp<-qlogis(pmin(pmax(p,1e-6),1-1e-6))
  a<-glm(y~1+offset(lp),family=binomial())
  s<-glm(y~lp,family=binomial())
  c(Intercept=unname(coef(a)[1]),Intercept_SE=sqrt(vcov(a)[1,1]),Slope=unname(coef(s)[2]),
    Brier=mean((p-y)^2),Scaled_Brier=1-mean((p-y)^2)/(mean(y)*(1-mean(y))))
}
model_path<-'output/cp2_final_model.rds'
model_md5<-unname(tools::md5sum(model_path))
frozen<-readRDS(model_path)
stopifnot(inherits(frozen,'workflow'))
selection<-read.csv('output/cp2_final_model_selection.csv')
stopifnot(selection$Model=='Random Forest')
rec<-workflows::extract_recipe(frozen,estimated=TRUE)
schema<-summary(rec,original=TRUE)
inputs<-as.character(schema$variable[schema$role=='predictor'])
align<-function(d) {
  for(nm in setdiff(inputs,names(d))) d[[nm]]<-if(grepl('_count$',nm)) 0 else NA_real_
  d[,c('RecordID',inputs),drop=FALSE]
}
predict_risk<-function(d) as.numeric(predict(frozen,new_data=align(d),type='prob')$.pred_Died)
message('Building the three sets using the shared Checkpoint 2 feature function')
sets<-setNames(lapply(c('a','b','c'),cp_build_table),c('a','b','c'))
ids<-lapply(sets,function(z) z$predictors$RecordID)
stopifnot(!length(intersect(ids$a,ids$b)),!length(intersect(ids$a,ids$c)),!length(intersect(ids$b,ids$c)))
# Preserve Checkpoint 2 column order and check every derived feature is present.
stopifnot(all(inputs %in% names(sets$a$predictors)),all(inputs %in% names(sets$b$predictors)),all(inputs %in% names(sets$c$predictors)))
bench_a<-sets$a$benchmarks
score_models<-lapply(c('SAPS-I','SOFA'),function(nm){
  x<-bench_a[[nm]];x[x==-1|!is.finite(x)]<-NA_real_;med<-median(x,na.rm=TRUE);x[is.na(x)]<-med
  list(model=glm(death~score,data=data.frame(death=bench_a$death,score=x),family=binomial()),median=med)
})
names(score_models)<-c('SAPS-I','SOFA')
saveRDS(score_models,'output/cp3_set_a_score_mappings.rds')
test<-bind_rows(lapply(c('b','c'),function(s){
  d<-sets[[s]]$predictors;d$p<-predict_risk(d);bm<-sets[[s]]$benchmarks
  stopifnot(identical(d$RecordID,bm$RecordID))
  for(nm in names(score_models)) {
    x<-bm[[nm]];x[x==-1|!is.finite(x)]<-NA_real_;d[[paste0(nm,'_missing')]]<-is.na(x)
    x[is.na(x)]<-score_models[[nm]]$median
    d[[nm]]<-as.numeric(predict(score_models[[nm]]$model,newdata=data.frame(score=x),type='response'))
  }
  mutate(d,Set=paste0('set-',s),.before=1)
}))
stopifnot(nrow(test)==8000,all(is.finite(test$p)),all(test$p>=0&test$p<=1))
cohorts<-test |> group_by(Set) |> summarise(Stays=n(),Deaths=sum(death),Death_rate=mean(death),Mean_risk=mean(p),SAPS_missing=sum(`SAPS-I_missing`),SOFA_missing=sum(SOFA_missing),.groups='drop')
save_table(cohorts,'cohorts')
save_table(test |> select(Set,RecordID,death,p,`SAPS-I`,SOFA,Age,Gender,ICUType),'predictions')
message('Computing discrimination and calibration')
auc_draws<-list()
performance<-bind_rows(lapply(c('set-b','set-c'),function(s){
  d<-filter(test,Set==s)
  bind_rows(lapply(c('p','SAPS-I','SOFA'),function(m){
    draws<-boot_auc(d$death,d[[m]])
    if(m=='p') auc_draws[[s]]<<-draws
    tibble(Set=s,Model=if(m=='p') 'Random Forest' else m,AUROC=auc(d$death,d[[m]]),Lower=quantile(draws,.025),Upper=quantile(draws,.975))
  }))
}))
save_table(performance,'discrimination')
cv<-read.csv('output/cp2_model_comparison.csv')
comparison<-performance |> left_join(cv |> select(Model,Mean_AUROC,SD_AUROC),by='Model')
save_table(comparison,'cv_comparison')
delta<-auc_draws[['set-b']]-auc_draws[['set-c']]
delta_ci<-quantile(delta,c(.025,.975))
rf_perf<-filter(performance,Model=='Random Forest')
roc_df<-bind_rows(lapply(split(test,test$Set),function(d) {
  r<-yardstick::roc_curve(tibble(y=factor(d$death,levels=c(0,1)),p=d$p),y,p,event_level='second')
  mutate(r,Set=d$Set[1])
}))
roc_plot<-ggplot(roc_df,aes(1-specificity,sensitivity,color=Set))+geom_path(linewidth=.9)+geom_abline(linetype=2,color='grey60')+coord_equal()+theme_minimal()+labs(x='False-positive rate',y='Sensitivity',title='Frozen random forest: discrimination')
save_plot(roc_plot,'roc')
b<-filter(test,Set=='set-b');cset<-filter(test,Set=='set-c')
pr<-yardstick::pr_curve(tibble(y=factor(b$death,levels=c(0,1)),p=b$p),y,p,event_level='second')
pr_auc<-yardstick::pr_auc_vec(factor(b$death,levels=c(0,1)),b$p,event_level='second')
pr_plot<-ggplot(pr,aes(recall,precision))+geom_path(color='#226699',linewidth=1)+geom_hline(yintercept=mean(b$death),linetype=2,color='#bb4444')+theme_minimal()+labs(title='Set-b precision-recall curve',subtitle=paste('Dashed baseline = death rate',fmt(mean(b$death))),x='Recall (sensitivity)',y='Precision (positive predictive value)')
save_plot(pr_plot,'precision_recall')
calibration<-bind_rows(lapply(split(test,test$Set),function(d) tibble::as_tibble_row(c(Set=d$Set[1],cal(d$death,d$p))))) |> mutate(across(-Set,as.numeric))
save_table(calibration,'calibration')
smooth_one<-function(d) {
  grid<-seq(min(d$p),max(d$p),length.out=101)
  smooth_predict<-function(x) {
    fit<-suppressWarnings(loess(death~p,data=x,span=.75,degree=1,control=loess.control(surface='direct',statistics='approximate')))
    pmin(1,pmax(0,as.numeric(predict(fit,newdata=data.frame(p=grid)))))
  }
  point<-smooth_predict(d)
  draws<-replicate(200,{i<-sample.int(nrow(d),nrow(d),TRUE);smooth_predict(d[i,])})
  tibble(Set=d$Set[1],Risk=grid,Observed=point,Lower=apply(draws,1,quantile,.025,na.rm=TRUE),Upper=apply(draws,1,quantile,.975,na.rm=TRUE))
}
curve<-bind_rows(lapply(split(test,test$Set),smooth_one))
save_table(curve,'calibration_curve')
cal_plot<-ggplot(curve,aes(Risk,Observed))+geom_ribbon(aes(ymin=Lower,ymax=Upper),fill='#77aadd',alpha=.35)+geom_line(color='#225588',linewidth=1)+geom_abline(linetype=2)+geom_rug(data=test,aes(x=p),inherit.aes=FALSE,alpha=.05,sides='b')+facet_wrap(~Set)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_minimal()+labs(title='Flexible calibration with 95% pointwise bootstrap bands',subtitle='Local linear loess, span 0.75; rug shows predicted risk distribution',x='Predicted probability',y='Observed death probability')
save_plot(cal_plot,'calibration',11,5.5)
risk_regions<-test |> mutate(Region=cut(p,c(0,.1,.2,.4,1),include.lowest=TRUE)) |> group_by(Set,Region) |> summarise(N=n(),Predicted=mean(p),Observed=mean(death),Gap=Observed-Predicted,.groups='drop')
save_table(risk_regions,'risk_regions')
threshold<-.15
threshold_metrics<-function(d) {
  flag<-d$p>=threshold;tp<-sum(flag&d$death==1);fp<-sum(flag&d$death==0)
  se<-tp/sum(d$death);ppv<-tp/(tp+fp)
  tibble(Set=d$Set[1],Threshold=threshold,TP=tp,FP=fp,FN=sum(!flag&d$death==1),TN=sum(!flag&d$death==0),Sensitivity=se,PPV=ppv,Flagged_per_100=100*mean(flag),Event1=min(se,ppv))
}
threshold_table<-bind_rows(lapply(split(test,test$Set),threshold_metrics));save_table(threshold_table,'threshold')
nb<-function(y,p,t) mean(p>=t & y==1)-mean(p>=t & y==0)*t/(1-t)
dca<-bind_rows(lapply(split(test,test$Set),function(d) bind_rows(lapply(seq(.05,.5,.01),function(t){
  tibble(Set=d$Set[1],Threshold=t,Strategy=c('Random Forest','SOFA','Treat all','Treat none'),Net_benefit=c(nb(d$death,d$p,t),nb(d$death,d$SOFA,t),nb(d$death,rep(1,nrow(d)),t),0))
}))))
save_table(dca,'decision_curve')
dca_plot<-ggplot(dca,aes(Threshold,Net_benefit,color=Strategy))+geom_line(linewidth=.9)+facet_wrap(~Set)+coord_cartesian(ylim=c(-.04,.15))+theme_minimal()+theme(legend.position='bottom')+labs(title='Decision curves: review at the end of the 48-hour window',subtitle='Treat-all values below -0.04 are clipped for readability',x='Risk threshold',y='Net benefit')
save_plot(dca_plot,'decision_curve',11,6)
nb15<-dca |> filter(abs(Threshold-threshold)<1e-8) |> select(-Threshold) |> pivot_wider(names_from=Strategy,values_from=Net_benefit)
save_table(nb15,'net_benefit_15')
b<-b |> mutate(Age_group=cut(Age,c(-Inf,50,65,80,Inf),labels=c('<=50','51-65','66-80','>80')),Sex=case_when(Gender=='0'~'Female',Gender=='1'~'Male',TRUE~'Missing'),ICU=case_when(ICUType=='1'~'Coronary care',ICUType=='2'~'Cardiac surgery recovery',ICUType=='3'~'Medical ICU',ICUType=='4'~'Surgical ICU',TRUE~'Missing'))
subgroups<-bind_rows(lapply(c('Age_group','Sex','ICU'),function(v) {
  g<-as.character(b[[v]]);g[is.na(g)]<-'Missing'
  bind_rows(lapply(unique(g),function(s){
    d<-b[g==s,];ok<-sum(d$death)>=10&&sum(d$death==0)>=10
    if(ok) {ci<-quantile(boot_auc(d$death,d$p,1000),c(.025,.975));cc<-cal(d$death,d$p)} else {ci<-c(NA,NA);cc<-rep(NA_real_,5);names(cc)<-names(cal(b$death,b$p))}
    tibble(Variable=v,Group=s,N=nrow(d),Deaths=sum(d$death),AUROC=if(ok) auc(d$death,d$p) else NA_real_,Lower=ci[1],Upper=ci[2],CITL=cc['Intercept'],CITL_lower=cc['Intercept']-1.96*cc['Intercept_SE'],CITL_upper=cc['Intercept']+1.96*cc['Intercept_SE'])
  }))
}))
save_table(subgroups,'subgroups')
sg_plot<-ggplot(subgroups,aes(AUROC,Group))+geom_errorbar(aes(xmin=Lower,xmax=Upper),orientation='y',width=.2)+geom_point(color='#226699')+facet_wrap(~Variable,scales='free_y',ncol=1)+theme_minimal()+labs(title='Set-b subgroup AUROC with bootstrap 95% intervals',x='AUROC',y=NULL)
save_plot(sg_plot,'subgroups',9,9)
message('Computing frozen-model explanations')
engine<-workflows::extract_fit_engine(frozen)
importance<-tibble(Predictor=names(engine$variable.importance),Importance=as.numeric(engine$variable.importance)) |> arrange(desc(Importance))
stopifnot(nrow(importance)>0)
save_table(importance,'importance')
imp_plot<-importance |> slice_head(n=15) |> ggplot(aes(Importance,reorder(Predictor,Importance)))+geom_col(fill='#226699')+theme_minimal()+labs(title='Set-a permutation importance in the frozen forest',y=NULL,x='Increase in prediction error after permutation')
save_plot(imp_plot,'importance')
# Age is chosen before examining test outcomes and has a matching CP2 GAM curve.
set.seed(3707)
bg<-sets$a$predictors[sample.int(4000,400),]
ages<-seq(20,90,2)
pdp<-tibble(Age=ages,Risk=vapply(ages,function(age){d<-bg;d$Age<-age;mean(predict_risk(d))},numeric(1)))
save_table(pdp,'age_pdp')
pdp_plot<-ggplot(pdp,aes(Age,Risk))+geom_line(color='#226699',linewidth=1)+geom_rug(data=bg,aes(Age),inherit.aes=FALSE,alpha=.15)+theme_minimal()+labs(title='Random forest: age partial dependence',subtitle='Average prediction across 400 fixed set-a background stays',y='Average predicted mortality risk')
save_plot(pdp_plot,'age_pdp')
cases<-bind_rows(b |> filter(death==0) |> arrange(desc(p),RecordID) |> slice_head(n=5) |> mutate(Case='Highest-risk survivor'),b |> filter(death==1) |> arrange(p,RecordID) |> slice_head(n=5) |> mutate(Case='Lowest-risk death'))
# Engine-space donor replacement is a sensitivity explanation, not SHAP.
bake_x<-function(d) recipes::bake(rec,new_data=align(d),all_predictors())
xa<-as.data.frame(bake_x(sets$a$predictors));xc<-as.data.frame(bake_x(cases))
ep<-function(x) as.numeric(predict(engine,data=x,num.threads=2)$predictions[,'Died'])
stopifnot(max(abs(ep(xc)-cases$p))<1e-10)
features<-head(importance$Predictor,30)
set.seed(4707);donors<-sample.int(nrow(xa),50)
local<-bind_rows(lapply(seq_len(nrow(cases)),function(i){
  bind_rows(lapply(features,function(f){
    x<-xc[rep(i,length(donors)),,drop=FALSE];x[[f]]<-xa[[f]][donors]
    tibble(RecordID=cases$RecordID[i],Predictor=f,Value=xc[[f]][i],Contribution=cases$p[i]-mean(ep(x)))
  }))
}))
save_table(local,'local_sensitivity')
driver_text<-function(id,high){
  z<-filter(local,RecordID==id) |> arrange(if(high) desc(Contribution) else Contribution) |> slice_head(n=3)
  paste(sprintf('%s = %.2f (%+.3f)',z$Predictor,z$Value,z$Contribution),collapse='; ')
}
raw_features<-sets$b$predictors[,inputs,drop=FALSE]
missing<-rowMeans(is.na(raw_features))
case_report<-cases |> transmute(RecordID,Case,Observed=death,Risk=p,Age,Missing_fraction=missing[match(RecordID,sets$b$predictors$RecordID)])
case_report$Drivers<-vapply(seq_len(nrow(cases)),function(i) driver_text(cases$RecordID[i],cases$death[i]==0),character(1))
case_report$Data_flags<-vapply(seq_len(nrow(cases)),function(i){
  d<-cases[i,];flags<-character()
  if(is.finite(d$Temp_min)&&d$Temp_min<30) flags<-c(flags,'Temp_min <30')
  if(is.finite(d$Temp_max)&&d$Temp_max>45) flags<-c(flags,'Temp_max >45')
  if(is.finite(d$HR_min)&&d$HR_min<=0) flags<-c(flags,'HR_min <=0')
  if(is.finite(d$RespRate_min)&&d$RespRate_min<=0) flags<-c(flags,'RespRate_min <=0')
  if(is.finite(d$GCS_min)&&(d$GCS_min<3||d$GCS_max>15)) flags<-c(flags,'GCS outside 3-15')
  if(is.finite(d$pH_min)&&(d$pH_min<6.5||d$pH_max>8)) flags<-c(flags,'pH outside 6.5-8')
  if(case_report$Missing_fraction[i]>.5) flags<-c(flags,'>50% predictor summaries missing')
  if(length(flags)) paste(flags,collapse='; ') else 'No selected range flag'
},character(1))
case_report$Judgment<-vapply(seq_len(nrow(cases)),function(i){
  lead<-filter(local,RecordID==cases$RecordID[i]) |> arrange(if(cases$death[i]==0) desc(Contribution) else Contribution) |> slice(1)
  context<-if(grepl('GCS',lead$Predictor)) 'The GCS contribution could reflect neurological status or unrecorded sedation context.' else if(grepl('Urine',lead$Predictor)) 'The urine-output contribution may reflect renal status, fluid management or collection practices that are not distinguished here.' else if(grepl('BUN',lead$Predictor)) 'The BUN contribution does not encode the underlying diagnosis, treatment response or subsequent renal trajectory.' else 'The dominant physiological summary omits its diagnosis and treatment context.'
  judgment<-if(case_report$Data_flags[i]!='No selected range flag') 'Data-quality/missingness contribution is plausible; flags require source review. Model error and unrecorded information remain possible.' else if(cases$death[i]==0) 'Model overestimates this observed outcome; recovery after severe early physiology is plausible, but cannot be confirmed.' else 'Model underestimates this observed outcome; later deterioration or omitted diagnosis/care limitations are plausible, but cannot be confirmed.'
  paste(judgment,context)
},character(1))
save_table(case_report,'confident_errors')
top_local<-local |> group_by(RecordID) |> slice_max(abs(Contribution),n=5,with_ties=FALSE) |> ungroup()
local_plot<-ggplot(top_local,aes(Contribution,reorder(Predictor,Contribution),fill=Contribution>0))+geom_col()+facet_wrap(~RecordID,ncol=2,scales='free_y')+scale_fill_manual(values=c('#4477aa','#cc6677'),guide='none')+theme_minimal()+labs(title='Local donor-replacement sensitivity for ten set-b errors',x='Original risk minus donor-replaced risk',y=NULL)
save_plot(local_plot,'local_errors',12,14)
stopifnot(identical(model_md5,unname(tools::md5sum(model_path))))
runtime<-as.numeric(difftime(Sys.time(),started,units='mins'))
save_table(data.frame(Minutes=runtime,Model_MD5=model_md5,R=R.version.string),'runtime')
capture.output(sessionInfo(),file='output/cp3_session_info.txt')
message('Checkpoint 3 analysis finished in ',round(runtime,2),' minutes')
