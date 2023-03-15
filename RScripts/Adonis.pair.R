adonis.pair<-function(dist.mat,Factor,nper=1000,corr.method="fdr"){
  require(vegan)
  comb.fact<-combn(levels(as.factor(Factor)),2)
  pv<-NULL
  R2<-NULL
  SS<-NULL
  MeanSqs<-NULL
  F.Model<-NULL
  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis2(as.dist(as.matrix(dist.mat)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]]) ~ Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]], permutations=nper)
    pv<-c(pv,model.temp$aov.tab[[6]][1])
    R2<-c(R2,model.temp$aov.tab$R2[1])
    SS<-c(SS,model.temp$aov.tab[[2]][1])
    MeanSqs<-c(MeanSqs,model.temp$aov.tab[[3]][1])
    F.Model<-c(F.Model,model.temp$aov.tab[[4]][1])
    }
  pv.corr<-p.adjust(pv,method=corr.method)
  data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),SumsOfSqs=SS,MeanSqs=MeanSqs,F.Model=F.Model,R2=R2,P.value=pv,P.value.corrected=pv.corr)
  }


adonis.pair<-function(dist.mat,Factor,nper=1000,corr.method="fdr"){
  require(vegan)
  comb.fact<-combn(levels(as.factor(metadata.nocontrols$plant.water)),2)
  pv<-NULL
  R2<-NULL
  SS<-NULL
  MeanSqs<-NULL
  F.Model<-NULL
  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis2(as.dist(as.matrix(Dist.clr)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]]) ~ Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]], permutations=nper)
    pv<-c(pv,model.temp$aov.tab[[6]][1])
    R2<-c(R2,model.temp$aov.tab$R2[1])
    SS<-c(SS,model.temp$aov.tab[[2]][1])
    MeanSqs<-c(MeanSqs,model.temp$aov.tab[[3]][1])
    F.Model<-c(F.Model,model.temp$aov.tab[[4]][1])
    }
  pv.corr<-p.adjust(pv,method=corr.method)
  data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),SumsOfSqs=SS,MeanSqs=MeanSqs,F.Model=F.Model,R2=R2,P.value=pv,P.value.corrected=pv.corr)
  }
