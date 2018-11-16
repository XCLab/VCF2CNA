trimed.median<-function(xx) {
  cutt<-0
  dnn<-quantile(xx,0.01)
  upp<-quantile(xx,0.99)
  if(sum(xx>dnn & xx<upp)>2) {
    use<-(xx>dnn & xx<upp)
    cutt<-quantile(xx[use],0.98)
  }
  return(as.numeric(cutt))
}
beutify.it.0<-function(xx,offsett=1e6) {
  zz<-ceiling((10*offsett+xx)/10^floor(log10(xx)))
  return(zz*(10^floor(log10(xx))))
}

plot.one.panel<-function(poss,MAF,CVG,medCVG,cnv,loh,xoffset,chrlen, nroww, upoffset.low,upoffset.high,nPrev,spacer,shrinkage=1,seq.type="WGS") {
  intervals<-floor(chrlen/1e7)
  MAF[MAF>=1]<-1
  MAF[MAF<=0]<-0
  if(!is.null(CVG)) {
    CVG[CVG>medCVG]<-medCVG
    if(length(poss)>0) { 
      for(j in 1:length(poss)) {
        lines(c(xoffset+poss[j],xoffset+poss[j]),c(upoffset.high,upoffset.high+shrinkage*CVG[j]/medCVG),col="orange")
      }
    }
  }
  cex.points<-0.4
  if(seq.type=="WES") cex.points<-1
  points(xoffset+poss,upoffset.high+shrinkage*MAF,pch=16,cex=cex.points,col="#555555")
  the.col<-"skyblue"
  if(nPrev%%2==1) the.col<-"springgreen"  
  rect(xleft=xoffset,xright=xoffset+chrlen,ytop=upoffset.high,ybottom=upoffset.low,col=the.col,border=NA)
  if(nPrev>0) {
    rect(xleft=xoffset-spacer*0.7,xright=xoffset-spacer*0.3,ytop=upoffset.high+shrinkage,ybottom=upoffset.low,col="red",border=NA)
  }
  for(j in 0:intervals) {
    lines(c(xoffset+j*1e7,xoffset+j*1e7),c(upoffset.low,upoffset.high),lwd=4)
    tick.len<-upoffset.high-upoffset.low
    for(k in 1:9) {
      if((j+k/10)*1e7<chrlen) {
        tick.height<-0.8
        if(k%%2==0) {
          tick.height<-0.5
        }
        lines(c(xoffset+(j+k/10)*1e7,xoffset+(j+k/10)*1e7),c(upoffset.low+tick.len*tick.height,upoffset.high),lwd=4)
      }
    }
  }
  lines(c(xoffset+min(poss),xoffset+max(poss)),c(upoffset.high+shrinkage*0.9,upoffset.high+shrinkage*0.9),lwd=1,col="gray",lty=2)
  lines(c(xoffset+min(poss),xoffset+max(poss)),c(upoffset.high+shrinkage*0.8,upoffset.high+shrinkage*0.8),lwd=3,col="gray")
  if(!is.null(cnv)) {
    for(j in 1:nrow(cnv)) {
      seg.mean<-cnv[j,"seg.mean"]
      if(seg.mean>0.25) {
        lines(c(xoffset+cnv[j,"loc.start"],xoffset+cnv[j,"loc.end"]),c(upoffset.high+shrinkage*0.95,upoffset.high+shrinkage*0.95),lwd=15,col="red")
      }
      if(seg.mean<0-0.25) {
        lines(c(xoffset+cnv[j,"loc.start"],xoffset+cnv[j,"loc.end"]),c(upoffset.high+shrinkage*0.85,upoffset.high+shrinkage*0.85),lwd=15,col="blue")
      }
    }
  }
  if(!is.null(loh)) {
    for(j in 1:nrow(loh)) {
      seg.mean<-loh[j,"seg.mean"]
      if(seg.mean>0.1) {
        lines(c(xoffset+loh[j,"loc.start"],xoffset+loh[j,"loc.end"]),c(upoffset.high+shrinkage*0.75,upoffset.high+shrinkage*0.75),lwd=15,col="purple")
      }
    }
  }
  xpos<-xoffset+(1:intervals)*1e7
  return(xpos)
}

arrange.panels<-function(chr.size.ifn,dt0,cnv,loh,maxcvg) {
  chr.sz<-read.table(chr.size.ifn,header=T,sep="\t",colClass=c("character",rep("numeric",3)))
  nSNP<-nrow(dt0)
  seq.type<-"WGS"
  if(nSNP<6e5) {
    seq.type<-"WES"
  }
  chrs<-apply(cbind("chr",c(1:22,"X","Y")),1,paste,collapse="")
  rows<-sort(unique(chr.sz[,"Row"]))
  max.len<-0
  panel.height<-1
  panel.dist<-0.1
  letter.track.height<-0.8
  track.height<-3*(panel.height+panel.dist)+letter.track.height
  spacer<-3e6
  shrinkage<-0.9
  for(i in 1:length(rows)) {
    len<-sum(chr.sz[chr.sz[,"Row"]==rows[i],"Length"])
    if(len>max.len) max.len<-len
  }  
  nroww<-length(rows)*track.height
  maxY<-nroww
  maxX<-max.len*1.05
  medd.D<-trimed.median(dt0[,"TumorTotal"])*1.5
  medd.G<-trimed.median(dt0[,"NormalTotal"])*1.5
  if(maxcvg!=-1) {
    medd.D<-maxcvg
    medd.G<-maxcvg
  }
  if(medd.D<1) medd.D<-1
  if(medd.G<1) medd.G<-1
  medd<-beutify.it.0(max(medd.D,medd.G),offsett=0)
  plot(-5,xlim=c(0,maxX),ylim=c(0,maxY),xlab="",ylab="",main="",xaxt="n",yaxt="n")
  for(i in 1:length(chrs)) {
    the.row<-chr.sz[chr.sz[,1]==chrs[i],"Row"]
    the.Y<-(length(rows)-the.row)*track.height
    starts<-0
    nPrev<-0
    nFut<-0
    if(i>1) {
      prev_chrs<-chrs[1:(i-1)]
      prev_chrs<-intersect(prev_chrs, chr.sz[chr.sz[,"Row"]==chr.sz[i,"Row"],"Chr"])
      nPrev<-length(prev_chrs)
      starts<-sum(chr.sz[as.character(chr.sz[,1]) %in% prev_chrs,"Length"])
      if(length(prev_chrs)>0)   starts<-starts+length(prev_chrs)*spacer
      fut_chrs<-setdiff(chrs[i:length(chrs)],chrs[i])
      fut_chrs<-intersect(fut_chrs, chr.sz[chr.sz[,"Row"]==chr.sz[i,"Row"],"Chr"])
      nFut<-length(fut_chrs)
      if(nFut==0) {
        themy.x<-starts+chr.sz[chr.sz[,"Chr"]==chrs[i],"Length"] + spacer
        panel.lbl<-c("D-G","G","D")
        yy<-(length(rows)-the.row)*track.height+letter.track.height+(0:2)*(panel.height+panel.dist)+panel.dist
        for(kk in 1:length(yy)) {
          lines(c(themy.x, themy.x), c(yy[kk],yy[kk]+shrinkage),lwd=3)
          for(ll in 0:5) {
            lines(c(themy.x,themy.x+maxX*0.002),c(yy[kk]+shrinkage*ll/5,yy[kk]+shrinkage*ll/5),lwd=3)
            text(themy.x+maxX*0.002,yy[kk]+shrinkage*ll/5,medd*ll/5,cex=3,pos=4)
          }
          text(themy.x+maxX*0.015,yy[kk]+shrinkage*2.5/5,panel.lbl[kk],cex=4,pos=4)
        }
      }
    }
    poss0<-dt0[as.character(dt0[,1])==chrs[i],"Pos"]
    ord<-order(as.numeric(poss0))
    poss<-poss0[ord]
    dCVG<-dt0[as.character(dt0[,1])==chrs[i],"TumorTotal"][ord]
    dMAF<-dt0[as.character(dt0[,1])==chrs[i],"dMAF"][ord]
    nMAF<-dt0[as.character(dt0[,1])==chrs[i],"nMAF"][ord]
    nCVG<-dt0[as.character(dt0[,1])==chrs[i],"NormalTotal"][ord]
    tmpcnv<-NULL
    tmploh<-NULL
    if(!is.null(cnv)) {
      if(chrs[i] %in% as.character(cnv$chrom)) {
        tmpcnv<-cnv[as.character(cnv$chrom)==chrs[i],]
      }
    }
    if(!is.null(loh)) {
      if(chrs[i] %in% as.character(loh$chrom)) {
        tmploh<-loh[as.character(loh$chrom)==chrs[i],]
      }
    }
    plot.one.panel(poss=poss, MAF=dMAF,CVG=dCVG,medCVG=medd,cnv=NULL,loh=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+2*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+2*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type)
    plot.one.panel(poss=poss, MAF=nMAF,CVG=nCVG,medCVG=medd,cnv=NULL,loh=NULL,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+1*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+1*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type)
    xpos<-plot.one.panel(poss=poss, MAF=abs(nMAF-dMAF),CVG=NULL,medCVG=medd,cnv=tmpcnv,loh=tmploh,xoffset=starts,
                   chrlen=chr.sz[chr.sz[,1]==chrs[i],"Length"], 
                   nroww=the.row, upoffset.low=the.Y+(letter.track.height+0*(panel.height+panel.dist)),
                   upoffset.high=the.Y+(letter.track.height+0*(panel.height+panel.dist))+panel.dist,
                   nPrev=nPrev,spacer=spacer,shrinkage=shrinkage,seq.type=seq.type)
    for(j in 1:length(xpos)) {
      text(xpos[j], the.Y+letter.track.height-0.1, labels=j*10,cex=3)
    }
    mid<-median(xpos)
    text(mid, the.Y+letter.track.height*0.5,chrs[i],cex=7)
  }
  for(j in 1:length(rows)) {
    the.x<-0-5e6
    yy<-(j-1)*track.height+letter.track.height+(0:2)*(panel.height+panel.dist)+panel.dist
    for(k in 1:length(yy)) {
      lines(c(the.x, the.x), c(yy[k],yy[k]+shrinkage),lwd=3)
      for(ll in 0:5) {
        lines(c(the.x-maxX*0.002,the.x),c(yy[k]+shrinkage*ll/5,yy[k]+shrinkage*ll/5),lwd=3)
        text(the.x-maxX*0.002,yy[k]+shrinkage*ll/5,ll/5,cex=3,pos=2)
      }
    }
  }
}

read.cnv<-function(cnv.ifn) {
  cnv<-read.table(cnv.ifn,header=T,colClass="character")
  thechr<-paste("chr",as.character(cnv[,1]),sep="")
  thechr<-gsub("chrchr","chr",thechr)
  thechr[thechr=="chr23" | thechr=="chrchr23"]<-"chrX"
  thechr[thechr=="chr24" | thechr=="chrchr24"]<-"chrY"
  ret<-data.frame(chrom=thechr,
                  loc.start=as.numeric(cnv[,"loc.start"]),
                  loc.end=as.numeric(cnv[,"loc.end"]),
                  num.mark=as.numeric(cnv[,"num.mark"]),
                  seg.mean=as.numeric(cnv[,"seg.mean"]),
                  LogRatio=as.numeric(cnv[,"LogRatio"]),
                  GMean=as.numeric(cnv[,"GMean"]),
                  DMean=as.numeric(cnv[,"DMean"]))
  return(ret)
}
read.loh<-function(loh.ifn) {
  loh<-read.table(loh.ifn,header=T,colClass="character")
  thechr<-paste("chr",as.character(loh[,1]),sep="")
  thechr<-gsub("chrchr","chr",thechr)
  thechr[thechr=="chr23" | thechr=="chrchr23"]<-"chrX"
  thechr[thechr=="chr24" | thechr=="chrchr24"]<-"chrY"
  ret<-data.frame(chrom=thechr,
                  loc.start=as.numeric(loh[,"loc.start"]),
                  loc.end=as.numeric(loh[,"loc.end"]),
                  num.mark=as.numeric(loh[,"num.mark"]),
                  seg.mean=as.numeric(loh[,"seg.mean"]))
  return(ret)
}


do.one.file<-function(chr.inf.ifn, dt.ifn, cnv.ifn,loh.ifn, jpg.ofn,maxcvg) {
  cnv<-read.cnv(cnv.ifn)
  loh<-read.loh(loh.ifn)
  dt0<-read.table(dt.ifn,header=T)
  if(!("dMAF" %in% colnames(dt0))) {
    dMAF<-dt0[,"TumorMutant"]/dt0[,"TumorTotal"]
    dt0<-cbind(dt0,dMAF)
  }
  if(!("nMAF" %in% colnames(dt0))) {
    nMAF<-dt0[,"NormalMutant"]/dt0[,"NormalTotal"]
    dt0<-cbind(dt0,nMAF)
  }
  xx<-dt0$NormalTotal
  if(maxcvg==-1) {
    if(sum(xx>15)>100) {
      upp<-as.numeric(quantile(xx[xx>10],0.98))
      cutofff<-mean(xx[xx<upp & xx>10]) + 4*sd(xx[xx<upp & xx>10])  
      dt0<-dt0[dt0$NormalTotal>15 & dt0$NormalTotal<cutofff,]
    }
    if(sum(xx>1)<5) {
      dt0<-dt0[dt0$TumorTotal>40,]
    }
  }
  if(maxcvg>1) {
    dt0<-dt0[dt0$NormalTotal>10,]
  }
#  jpeg(jpg.ofn,width=1200,height=6000)
  jpeg(jpg.ofn,width=8000,height=4000)
  arrange.panels(chr.size.ifn=chr.inf.ifn,dt0,cnv,loh,maxcvg)
  dev.off()
}


args<-commandArgs(TRUE)
if(length(args)==5) {
  do.one.file(chr.inf.ifn=args[1], dt.ifn=args[2], cnv.ifn=args[3], loh.ifn=args[4], jpg.ofn=args[5],maxcvg=-1)
}

if(length(args)==6) {
  do.one.file(chr.inf.ifn=args[1], dt.ifn=args[2], cnv.ifn=args[3], loh.ifn=args[4], jpg.ofn=args[5],maxcvg=as.integer(args[6]))
}

if(length(args)!=5 & length(args)!=6) {
  print("Usage: Rscript ~ <chr.size.ifn> <data.ifn> <cnv.ifn> <loh.ifn> <jpg.ofn> | <max_cvg>")
}

