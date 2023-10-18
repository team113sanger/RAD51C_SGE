pdf(file='LIG4_vs_WT.pdf',width=4.5,height=4.5);
gstable=read.table('LIG4_vs_WT.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ARRDC3","GSK3B","LOC100289561","TAOK1","TMA7","hsa-mir-6727","ERH","hsa-mir-483","MRPL40","PDE4DIP")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='LIG4_rep1,LIG4_rep2,LIG4_rep3_vs_WT_rep1,WT_rep2,WT_rep3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3.82122109335454,6.07428068332832,20.8262505674364,1.07233648498053,3.33316950066674,2.9298914381715),c(16.072918295394,23.308286343004,35.4992907399485,11.1925120619843,4.91305331023132,8.70655682690679),c(36.1752175880607,25.9922708309863,48.7681367454137,24.2951234878401,16.1263749829946,8.99746803353375),c(36.1900861526263,32.507996331417,37.2505890831193,11.7286803044745,13.8721505230061,25.7040830426819),c(16.3702895867056,18.8232070012442,35.6728428280105,19.3020567296495,3.56437200938351,4.82081428124672),c(2.79529013832939,0.706311707363759,1.92485043123276,0.234573606089491,2.29275821144128,3.96886003326777))
targetgene="ARRDC3"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5.67979166405228,5.38562676864866,3.92858817522096,6.23295581894932,2.96709886186519,3.30392013240616),c(13.6642108357697,11.9190100617634,22.3408869723409,5.328171909747,9.78757286900985,8.20785190126058),c(30.0939746807377,18.3111310134054,26.5692469360325,9.08134960717885,6.70487275278627,4.92471114075635),c(34.7180982606337,24.7915409284679,33.2115586700407,7.07071869784036,5.0864551917689,10.2857890914531),c(27.1351303321869,25.5155104285158,19.1222846119189,49.6123176879273,14.5272242977036,18.140391670381),c(14.4819818868767,12.5193750130226,15.3041386745556,4.8422694399902,9.07469846713315,2.36884839681951))
targetgene="GSK3B"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(14.794221742754,14.5853367570616,13.6632825692424,14.0409058502138,8.05355405363409,13.1325630420169),c(0.95158813219724,0.476760402470537,1.00975760326965,0.0167552575778208,0.0385337514527947,0.166234975215404),c(6.51243127972486,8.38745152494463,8.48827485248546,1.39068637895912,0.732141277603099,5.17406360357945),c(1.33817081090237,2.24253967087993,1.79863073082406,0.100531545466925,0.211935632990371,0.644160528959691),c(0.892113873934913,0.123604548788658,0.378659101226117,0.0502657727334623,0.0385337514527947,1.62079100835019),c(2.48305028245217,1.92469940256624,2.09840251929473,0.251328863667311,0.250469384443165,4.32210935560051))
targetgene="LOC100289561"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(12.2070915083427,8.4757404883651,14.0892740581218,15.3143054261282,15.3749668296651,7.62602948800667),c(11.6123489257195,14.0026295984865,19.7376056514114,2.61382018214004,2.8900313589596,11.6364482650783),c(22.4366639294631,28.4290462213913,48.0265960055125,12.0135196832975,11.9454629503664,13.8390616866824),c(5.5757117120932,2.89587800019141,7.69940172493105,0.552923500068085,0.481671893159933,1.39221791742901),c(9.97680682350544,16.5630095376801,14.0577191330196,3.552114606498,2.73589635314842,5.15328423167753),c(11.4339261509325,11.0714360129269,5.26967249206346,6.88641086448433,11.1940547970369,11.4702132898629))
targetgene="TAOK1"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.327108420442801,1.05946756104564,0.978202678167469,0.0837762878891038,0.211935632990371,0.041558743803851),c(0.312239855877219,0.529733780522819,0.583766114390264,0.0167552575778208,0.0192668757263973,0.0207793719019255),c(5.05531195229784,3.61984750023926,5.91654845665808,0.569678757645906,1.21381317076303,2.36884839681951),c(1.2340908589433,2.17190850014356,1.65663356786426,0.201063090933849,0.327536887348755,2.18183404970218),c(0.341976985008383,0.26486689026141,1.16753222878053,0.0502657727334623,0.0770675029055894,22.2962660507661),c(1.1448794715498,0.988836390309262,0.88353790286094,0.184307833356028,1.05967816495185,0.0207793719019255))
targetgene="TMA7"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.208159903918146,0.423787024418255,0.205107013164147,0.0335105151556415,0.0770675029055894,0.0),c(1.62067353764842,2.17190850014356,5.36433726736999,0.150797318200387,0.134868130084781,0.685719272763542),c(4.25240946575642,2.41911759772087,4.13369518838511,0.385370924289878,0.269736260169563,3.96886003326777),c(3.91043248074803,1.7481214757253,1.53041386745556,0.234573606089491,0.154135005811179,0.664939900861617))
targetgene="hsa-mir-6727"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3.00345004224754,4.25552803686665,3.3448220608307,0.485902469756802,0.308270011622358,3.55327259522926),c(11.0027377785306,10.1885463787222,8.18850306401478,4.69147212178981,2.60102822306364,2.30651028111373),c(0.758296792844676,0.0353155853681879,0.962425215616381,0.0335105151556415,0.0963343786319867,0.145455603313479),c(3.95503817444478,5.13841767107134,3.13971504766655,0.686965560690651,0.905543159140675,1.37143854552708),c(0.23789703304931,0.229551304893222,0.536433726736999,0.0335105151556415,0.0963343786319867,0.207793719019255),c(0.252765597614892,0.0529733780522819,0.489101339083735,0.0502657727334623,0.0,0.0207793719019255))
targetgene="ERH"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.431188372401875,0.388471439050067,0.347104176123941,0.067021030311283,0.057800627179192,0.18701434711733),c(0.223028468483728,0.759285085416041,0.347104176123941,0.067021030311283,0.0385337514527947,0.0207793719019255),c(0.163554210221401,0.794600670784229,0.0473323876532646,0.0167552575778208,0.173401881537576,0.041558743803851),c(0.0594742582623275,0.353155853681879,0.236661938266323,0.0167552575778208,0.0385337514527947,0.0831174876077021))
targetgene="hsa-mir-483"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.95158813219724,0.0882889634204698,1.89329550613058,0.117286803044745,0.0192668757263973,0.145455603313479),c(2.0518619100503,2.18956629282765,4.02325295052749,0.251328863667311,0.57800627179192,2.57664211583876),c(23.3139092388324,21.5425070745946,24.3761796414313,9.60076259209129,51.9627638340936,13.7351648271728),c(1.65041066677959,1.7657792684094,0.946647753065293,0.418881439445519,0.366070638801549,3.05456766958305),c(0.282502726746056,0.26486689026141,0.599543576941352,0.0502657727334623,0.0,0.0),c(0.0446056936967456,0.229551304893222,0.678430889696793,0.0167552575778208,0.0,0.0831174876077021))
targetgene="MRPL40"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2.33436463679636,2.86056241482322,4.05480787562967,0.351860409134236,0.366070638801549,0.394808066136585),c(5.93255726166717,5.59752028085779,6.59497934635487,7.87497106157576,3.31390262494034,2.16105467780025),c(18.7641284817643,11.089093805611,14.4206007716946,3.23376471251941,5.7607958421928,9.08058552114145),c(0.862376744803749,1.13009873178201,2.09840251929473,0.067021030311283,0.0963343786319867,0.207793719019255),c(1.97751908722239,3.83174101244839,2.09840251929473,0.335105151556415,0.154135005811179,0.374028694234659),c(0.743428228279094,0.282524682945503,0.930870290514205,0.134042060622566,0.057800627179192,0.540263669450064))
targetgene="PDE4DIP"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("OSTN","RFPL3","ATN1","hsa-mir-5697","TCP10L2","TCERG1L","MMP11","L3MBTL2","ZIM3","MORC2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='LIG4_rep1,LIG4_rep2,LIG4_rep3_vs_WT_rep1,WT_rep2,WT_rep3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(17.530037622821,18.8232070012442,19.1065071493678,29.4557428218089,30.3645961448022,14.5663397032498),c(5.93255726166717,7.50456189073994,6.59497934635487,13.8733532744356,5.97273147518318,9.51695233108189),c(13.7088165294665,13.6671315374887,11.2966631865792,34.0969491708652,9.22883347294433,19.0962427778695),c(7.16664812061047,12.9078464520727,13.600172719038,7.69066322821973,18.1686638099927,19.6988445630254),c(5.81360874514251,7.68113981758088,3.42370937358614,18.2632307598246,6.26173461107914,8.8727918021222),c(14.258953418393,15.8566978303164,14.5152655470011,29.8746242612544,29.4783198613879,12.6338581163707))
targetgene="OSTN"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(13.7534222231632,14.4264166229048,11.3439955742324,10.187196607315,91.7295953333777,18.3689647613022),c(9.87272687154637,13.8083938789615,10.9811139355574,5.52923500068085,10.8472510339617,5.98445910775455),c(22.9124579955617,27.4225520383979,26.0801455969488,68.8305981296877,15.1822980724011,69.0290734581966),c(5.96229439079833,5.80941379306692,11.9277616886227,17.6600414870231,4.12311140544903,19.0131252902619),c(11.1514234241864,17.6931082694622,10.6340097594335,20.4414142449413,26.5497547509755,34.9093447952349),c(7.67217931584025,9.55286584209484,7.17874546074514,7.17125024330728,5.10572206749529,9.47539358727804))
targetgene="RFPL3"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(4.29701515945316,6.86888135411255,6.04276815706679,1.89334410629375,16.1071081072682,13.9637379180939),c(10.0808867754645,8.75826517131061,18.8067353608972,9.24890218295706,20.3265538913492,51.1588136225406),c(3.12239855877219,3.01948254898007,2.03529266909038,6.91992137963997,3.02489948904438,6.1922528267738),c(4.20780377205967,11.5835120007656,4.96990070359279,11.963253910564,13.3712117541198,4.1558743803851),c(30.510294488574,28.6409397336004,36.5879356559736,57.956435961682,43.0036666213189,26.0365529931127),c(6.94361965212674,5.98599171990786,7.11563561054078,7.53986591001934,10.9435854125937,4.61302056222747))
targetgene="ATN1"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(20.1022992926667,17.5871615133576,15.8405724012926,11.4605961832294,48.8800637178701,28.0937108114033),c(12.5639370579167,11.6364853788179,9.3718127553464,18.3134965325581,13.8914173987325,18.6390965960272),c(5.81360874514251,5.01481312228269,9.07204096687572,16.2358445929083,12.7161379794222,7.10654519045853),c(4.46056936967456,4.11426569539389,5.33278234226782,11.6449040165854,10.3655791408018,13.7351648271728))
targetgene="hsa-mir-5697"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6.58677410255277,3.86705659781658,7.51007217431799,1.02207071224707,1.61841756101738,4.73769679363902),c(4.31188372401874,4.41444817102349,8.04650590105499,13.3539402895231,3.31390262494034,10.4520240666685),c(16.9204264756322,18.0286063304599,16.81877507946,75.7002537365942,50.3058125216235,43.2003141841032),c(4.69846640272387,4.78526181738946,4.43346697685579,6.75236880386177,18.3227988158039,21.8598992408257),c(11.4785318446292,11.830721098343,12.8270770540347,5.26115087943572,9.63343786319867,14.7325746784652))
targetgene="TCP10L2"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(7.77625926779932,5.13841767107134,10.8233393100465,8.71273394046679,5.70299521501361,23.2313377863527),c(3.76174683509221,4.46742154907577,4.13369518838511,10.4552807285601,1.04041128922546,9.0390267773376),c(37.766153996578,49.0003746983608,46.5908469133635,51.4553960214875,90.1882452752659,72.7070222848374),c(12.266565766605,11.3363029031883,13.5055079437315,38.3360293380539,29.3819854827559,15.3975145793268),c(2.39383889505868,1.27136107325477,1.92485043123276,2.99919110642992,0.385337514527947,0.872733619880872),c(27.9082956895972,23.4495486844768,31.0342688379905,45.7753637026063,53.0802426262247,50.3484181183655))
targetgene="TCERG1L"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(10.5566808415631,11.7600899276066,11.7857645256629,11.778946077208,61.9815392118202,24.3326444971548),c(11.5677432320227,14.1968653180116,9.65580708126599,18.1124334416242,30.114126760359,46.878263010744),c(0.654216840885603,0.635680536627383,1.24641954153597,0.067021030311283,0.0192668757263973,2.70131834725032),c(10.6161550998255,8.08726904931504,4.5439092147134,28.5844694277622,37.8786776780972,10.5767002980801),c(2.45331315332101,2.63101110993,1.35686177939359,0.43563669702334,0.655073774697509,4.11431563658125),c(23.968126079718,22.4607122941675,28.5729846800207,8.34411827375474,22.2147077125361,25.5170686955645))
targetgene="MMP11"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(11.9989316044246,7.327983963899,15.3514710622088,11.9129881378306,7.64894966337974,15.1066033726998),c(12.9505197366218,12.0779301959203,10.555122446678,8.96406280413411,23.7367908949215,16.3741450587173),c(14.1251363373028,14.5500211716934,5.83766114390264,70.5899001753588,77.6647760531077,12.0728150750187),c(10.631023664391,12.3251392934976,8.44094246483219,31.1647790947466,7.30214590030459,14.7533540503671),c(24.607474356038,19.4588875378716,23.6346389015301,18.0956781840464,19.9219495010949,26.8677278691897),c(6.22992855297881,3.2313760611892,5.64833159328958,3.18349893978594,16.0878412315418,4.44678558701206))
targetgene="L3MBTL2"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(0.847508180238167,2.47209097577315,0.362881638675029,0.603189272801547,3.54510513365711,1.35065917362516),c(13.7385536585977,15.0091237814799,13.3319558556695,16.9730759263324,17.8796606740967,24.7482319351933),c(30.7630600861889,27.2106585261888,25.3228273944966,59.7492485225088,67.1450619064947,27.1378597039147),c(5.81360874514251,9.69412818356759,7.65206933727778,16.822278608132,24.4689321725246,34.4937573571964),c(10.6161550998255,4.9794975369145,9.84513663187905,3.51860409134236,3.96897639963785,2.90911206626957),c(5.18912903338807,8.31682035420826,7.79406650023758,3.48509357618672,5.72226209074001,9.24682049635686))
targetgene="ZIM3"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(8.02902486541421,8.97015868351973,7.47851724921581,5.59625603099213,2.79369698032761,12.4052850254495),c(10.3782580667761,9.80007493967215,11.5017701997433,18.8999305477818,11.6179260630176,20.1767701167697),c(16.7866093945419,18.2758154280373,13.8526121198555,14.6273398654375,19.3054094778501,16.9767468438732),c(12.7423598327037,12.395770464234,16.0614568770078,26.456551715379,27.4167641586634,33.7456999687271),c(3.2859527689936,6.30383198822155,2.87149818429805,16.3363761383752,2.21569070853569,14.0052966618978),c(16.2067353764843,24.049913635736,14.4363782342457,25.0993758515755,29.4783198613879,18.4728616208118))
targetgene="MORC2"
collabel=c("WT_rep1","WT_rep2","WT_rep3","LIG4_rep1","LIG4_rep2","LIG4_rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("LIG4_vs_WT_summary.Rnw");
library(tools);

texi2dvi("LIG4_vs_WT_summary.tex",pdf=TRUE);

