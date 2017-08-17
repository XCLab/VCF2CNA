require(tree)||stop("Package tree not available. Please install the package")
DCODE="D"
GCODE="G"
gap.norm=T
mapability = T
withXY = T
merge.diff= T
quality.merge = T
thresh.diff=0.125
GMean.thresh = 0.2
gap.thresh = 1000
LOH.chr = 1:22
forced=F
plotRaw = F
toPerm = F
toRedraw = F
toNorm = F
plotAI = T
g.setPath = F
dev = 1e-6
window=100
disease = "SJCBF"
setGC = F
setMap = F
chr.toScale = array(0, dim=0)
chr.scale = array(0, dim = 0)
drawDetail = F
setPThresh=F
toDraw=F
detail.chr = -1
detail.start= -1
detail.stop = -1
detail.ymax = 2
SV.toFilter = F
SV.filter = 3
setPThresh = F
reDoSV = F
hg18 = T
GRCh38 = F
noLogWeight = T
singleBAM = F
single="G"
debug=F
unfilteredSV=T
args=(commandArgs(TRUE))
if (length(args) > 0) for (i in 1:length(args)) eval(parse(text=args[[i]]))
ts = paste("_", format(Sys.time(), "%Y%m%d%H%M%S"), sep="")
logthresh.diff=log2(1+thresh.diff)
auto.ref = !(exists("norm.chr") | exists("norm.seg") | singleBAM)
print(paste("Auto.ref", auto.ref))
Male=F
if (GRCh38) {
	mapability = F
	gap.norm = T
}

fun.help = function(mark.tmp, p.thresh)
{
  t = try(tree(diff~ChrPos, data=mark.tmp, mindev=dev), silent=T)
  de = dev
  while (class(t) == "try-error")
  {
    de = de * 2
    t = try(tree(diff~ChrPos, data=mark.tmp, mindev=de), silent=T)
  }
  t.p = try(prune.tree(t), silent=T)
  if (class(t.p) == "try-error")
  {
    t.p.best = t
  }
  else
  {
    t.p.bic = length(mark.tmp$diff)*log(t.p$dev/length(mark.tmp$diff)) + log(length(mark.tmp$diff)) * (t.p$size + 1)
    while (which(t.p.bic==min(t.p.bic))[1] == 1 & de > dev/1000)
    {
      t1 = t
      t1.p = t.p
      de = de / 2
      t = try(tree(diff~ChrPos, data=mark.tmp, mindev=de), silent=T)
      if (class(t) == "try-error")
      {
        t = t1
        t.p = t1.p
        break
      }
      t.p = prune.tree(t)
      t.p.bic = length(mark.tmp$diff)*log(t.p$dev/length(mark.tmp$diff)) + log(length(mark.tmp$diff)) * (t.p$size + 1)
    }
	
    t.p.best = prune.tree(t, best= min(t.p$size[which(t.p.bic==min(t.p.bic, na.rm=T))]))
  }

  seg.chr = unique(t.p.best$where)
  for (zz in 1:length(seg.chr)) seg.chr[zz] = sum(t.p.best$where==seg.chr[zz])
    return(seg.chr)
}

get_threshold = function(diff1, log1, diff2, log2, threshold)
{
  final_threshold = threshold
  abs_diff = abs(diff1 - diff2)
  abs_log  = abs(log1 - log2)
  
  if( abs_diff > 8 & abs_log > 3)
  {
    final_threshold = threshold * 1e6
  }
  return(final_threshold)
}

find_merge_segments = function(mark.chr, step, overlap, p.thresh, type, debug)
{
   current = 1
   offset = 0
   offset2 = 0
   num.mark1 = -1
   while (current <= nrow(mark.chr))
   {
      if (current + step*1.5 - 1 < nrow(mark.chr))
      {
         mark.tmp = mark.chr[current + ((-offset):(step - 1)),]
         current = current + step
      } else {
         mark.tmp = mark.chr[(current - offset):nrow(mark.chr),]
         current = nrow(mark.chr) + 1
      }
      if( type == "log")
      {
        mark.tmp$diff = mark.tmp$Log
      }
      mark.seg = fun.help(mark.tmp, p.thresh=p.thresh)
      if (debug) print(paste("Finishing", type, current, "bins.", date()))
      if (num.mark1[1] < 0) {
	 num.mark1 = mark.seg
      } else {
	 mark.seg[1] = mark.seg[1] + offset2
	 if (length(num.mark1) > 1) {
	    num.mark1 = c(num.mark1[1:(length(num.mark1) - 1)], mark.seg)
	 } else {
	    num.mark1 = mark.seg
	 }
      } 
      if (num.mark1[length(num.mark1)] > overlap) {
	 offset = overlap
      } else {
	 offset = num.mark1[length(num.mark1)]
      }
      offset2 = num.mark1[length(num.mark1)] - offset
   }

   Start = array(1, dim=length(num.mark1))
   if (length(Start) > 1)
   {
      for (i in 2 : length(num.mark1)) Start[i] = Start[i - 1] + num.mark1[i - 1]
      p.val = array(1, length(num.mark1)-1)
      for (i in 1: length(p.val))
      {
         p.val[i] = t.test(mark.chr$diff[Start[i]:(Start[i] + num.mark1[i] - 1)], mark.chr$diff[Start[i + 1]:(Start[i + 1] + num.mark1[i + 1] - 1)])$p.value
      }
      t.lim = array(1, length(p.val))
      for ( i in 1: length(t.lim))
      {
         diff_sig1 = mean(mark.chr$diff[Start[i]:(Start[i] + num.mark1[i] - 1)])
         log_sig1  = mean(mark.chr$Log[Start[i]:(Start[i] + num.mark1[i] - 1)])
         diff_sig2 = mean(mark.chr$diff[Start[i + 1]:(Start[i + 1] + num.mark1[i + 1] -1)])
         log_sig2  = mean(mark.chr$Log[Start[i + 1]:(Start[i + 1] + num.mark1[i + 1] -1)])
         t.lim[i]  = get_threshold(diff_sig1, log_sig1, diff_sig2, log_sig2, p.thresh) 
      }
      p.max.id = which(p.val == max(p.val))[1]
      
      while (!all(p.val < t.lim))
      {
         if ( p.val[p.max.id] > t.lim[p.max.id])
         {
            num.mark1[p.max.id] = num.mark1[p.max.id] + num.mark1[p.max.id+1]
            num.mark1 = num.mark1[-(p.max.id + 1)]
            p.val = p.val[-p.max.id]
            t.lim = t.lim[-p.max.id]
            Start = Start[-(p.max.id+1)]
            if (p.max.id > 1)
            {
               p.val[p.max.id-1] = t.test(mark.chr$diff[Start[p.max.id-1]:(num.mark1[p.max.id - 1] + Start[p.max.id-1] - 1)], mark.chr$diff[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)])$p.value
               diff_sig1         = mean(mark.chr$diff[Start[p.max.id-1]:(num.mark1[p.max.id - 1] + Start[p.max.id-1] - 1)])
               log_sig1          = mean(mark.chr$Log[Start[p.max.id-1]:(num.mark1[p.max.id - 1] + Start[p.max.id-1] - 1)])
               diff_sig2         = mean(mark.chr$diff[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)])
               log_sig2          = mean(mark.chr$Log[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)])
               t.lim[p.max.id-1] = get_threshold(diff_sig1, log_sig1, diff_sig2, log_sig2, p.thresh) 
            }
            if (p.max.id <= length(p.val))
            {
               p.val[p.max.id] = t.test(mark.chr$diff[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)], mark.chr$diff[Start[p.max.id+1]:(num.mark1[p.max.id+1] + Start[p.max.id+1] - 1)])$p.value
               diff_sig1 = mean(mark.chr$diff[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)])
               log_sig1 = mean(mark.chr$Log[Start[p.max.id]:(num.mark1[p.max.id] + Start[p.max.id] - 1)])
               diff_sig2 = mean(mark.chr$diff[Start[p.max.id+1]:(num.mark1[p.max.id+1] + Start[p.max.id+1] - 1)])
               log_sig2 = mean(mark.chr$Log[Start[p.max.id+1]:(num.mark1[p.max.id+1] + Start[p.max.id+1] - 1)])
               t.lim[p.max.id] = get_threshold(diff_sig1, log_sig1, diff_sig2, log_sig2, p.thresh)
            }  
            p.max.id = which(p.val == max(p.val))[1]
         }
         else
         {
           p.val[p.max.id] = 0
           p.max.id = which(p.val == max(p.val))[1]
         }
      }   
   }
   return(num.mark1)
}  

fun.chr = function(chr, step=20000, overlap = round(step/20), p.thresh=1e-9, debug=F) {
        mark.chr = cn.D[which(cn.D$SoftChr == chr),]
	if (sum(cn.D$SoftChr == chr) == 0) return(NULL)
	chr2 = mark.chr$Chr[1]

        num.mark1 = find_merge_segments(mark.chr, step, overlap, p.thresh, "diff", debug)
        num.mark2 = find_merge_segments(mark.chr, step, overlap, p.thresh, "log", debug)

	n1 = 1
	n2 = 1
	n3 = 1
	num.mark = array(0, length(c(num.mark1, num.mark2)))
	l1 = num.mark1[1]
	l2 = num.mark2[1]
	while ((n1 <= length(num.mark1)) & (n2 <= length(num.mark2))) {
	    if (l1 <= l2) {
		num.mark[n3] = l1
		l2 = l2 - l1
		l1 = 0
	    } else {
		num.mark[n3] = l2
		l1 = l1 - l2
		l2 = 0
	    }
	    if (l1 == 0) {
		n1 = n1 + 1
		if (n1 <= length(num.mark1)) l1 = num.mark1[n1]
	    }
	    if (l2 == 0) {
		n2 = n2 + 1
		if (n2 <= length(num.mark2)) l2 = num.mark2[n2]
	    }
	    n3 = n3 + 1
	}
	if (n1 < length(num.mark1)) {
	    if (l1 > 0) {
		num.mark[n3] = l1
		n3 = n3 + 1
		n1 = n1 + 1
	    }
	    while (n1 <= length(num.mark1)) {
		num.mark[n3] = num.mark1[n1]
		n3 = n3 + 1
		n1 = n1 + 1
	    }
	} else if (n2 < length(num.mark2)) {
	    if (l2 > 0) {
		num.mark[n3] = l2
		n3 = n3 + 1
		n2 = n2 + 1
	    }
	    while (n2 <= length(num.mark2)) {
		num.mark[n3] = num.mark2[n2]
		n3 = n3 + 1
		n2 = n2 + 1
	    }
	} else {
	    if (l1 > 0) {
		num.mark[n3] = l1
		n3 = n3 + 1
	    }
	    if (l2 > 0) {
		num.mark[n3] = l2
		n3 = n3 + 1
	    }
	}
	num.mark = (num.mark[1:(n3-1)])
	seg.chr =  data.frame(matrix(0, nrow=length(num.mark), ncol=9, dimnames = list(NULL, c("chrom", "loc.start", "loc.end", "num.mark", "length.ratio", "seg.mean", "GMean", "DMean", "LogRatio"))))
	seg.chr$num.mark = num.mark
	seg.chr$chrom = chr2
	Start = array(1, dim=nrow(seg.chr))
	if (nrow(seg.chr) > 1) {
	    for (i in 2 : nrow(seg.chr)) Start[i] = Start[i - 1] + seg.chr$num.mark[i - 1]
	}
	seg.chr$loc.start = mark.chr$ChrPos[Start] - 49
	seg.chr$loc.end = mark.chr$ChrPos[Start + seg.chr$num.mark - 1] + 50
	seg.chr$GMean = Start
	seg.chr$DMean = Start
	for (j in 1:length(Start)) {
	    seg.chr$seg.mean[j] = mean(mark.chr$diff[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	    seg.chr$LogRatio[j] = mean(mark.chr$Log[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	    seg.chr$GMean[j] = mean(mark.chr$GMean[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	    seg.chr$DMean[j] = mean(mark.chr$Mean[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	}

	seg.diff = abs(seg.chr$seg.mean[-nrow(seg.chr)] - seg.chr$seg.mean[-1])
	seg.diff.id = which(seg.diff == min(seg.diff))[1]
	while (min(seg.diff) < thresh.diff) {
	    seg.chr$seg.mean[seg.diff.id] = (seg.chr$seg.mean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$seg.mean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$GMean[seg.diff.id] = (seg.chr$GMean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$GMean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$DMean[seg.diff.id] = (seg.chr$DMean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$DMean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$LogRatio[seg.diff.id] = (seg.chr$LogRatio[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$LogRatio[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$num.mark[seg.diff.id] = seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1]
	    seg.chr$loc.end[seg.diff.id] = seg.chr$loc.end[seg.diff.id + 1]
	    seg.chr = seg.chr[-(seg.diff.id + 1),]
	    seg.diff = seg.diff[-seg.diff.id]
	    if (seg.diff.id > 1)
		seg.diff[seg.diff.id-1] = abs(seg.chr$seg.mean[seg.diff.id] - seg.chr$seg.mean[seg.diff.id-1])
	    if (seg.diff.id <= length(seg.diff))
		seg.diff[seg.diff.id] = abs(seg.chr$seg.mean[seg.diff.id+1] - seg.chr$seg.mean[seg.diff.id])
	    seg.diff.id = which(seg.diff == min(seg.diff))[1]
	}

	seg.diff = abs(seg.chr$LogRatio[-nrow(seg.chr)] - seg.chr$LogRatio[-1])
	seg.diff.id = which(seg.diff == min(seg.diff))[1]
	while (min(seg.diff) < logthresh.diff) {
	    seg.chr$seg.mean[seg.diff.id] = (seg.chr$seg.mean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$seg.mean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$GMean[seg.diff.id] = (seg.chr$GMean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$GMean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$DMean[seg.diff.id] = (seg.chr$DMean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$DMean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$LogRatio[seg.diff.id] = (seg.chr$LogRatio[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$LogRatio[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$num.mark[seg.diff.id] = seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1]
	    seg.chr$loc.end[seg.diff.id] = seg.chr$loc.end[seg.diff.id + 1]
	    seg.chr = seg.chr[-(seg.diff.id + 1),]
	    seg.diff = seg.diff[-seg.diff.id]
	    if (seg.diff.id > 1)
		seg.diff[seg.diff.id-1] = abs(seg.chr$LogRatio[seg.diff.id] - seg.chr$LogRatio[seg.diff.id-1])
	    if (seg.diff.id <= length(seg.diff))
		seg.diff[seg.diff.id] = abs(seg.chr$LogRatio[seg.diff.id+1] - seg.chr$LogRatio[seg.diff.id])
	    seg.diff.id = which(seg.diff == min(seg.diff))[1]
	}

	seg.chr$length.ratio = seg.chr$num.mark * 100 / (seg.chr$loc.end - seg.chr$loc.start + 1)
	return(seg.chr)
}

fun.chr.germline = function(chr, step=20000, overlap = round(step/20), p.thresh=1e-9, debug=F) {
	mark.chr = cn.D[which(cn.D$SoftChr == chr),]
	chr2 = mark.chr$Chr[1]
        num.mark1 = find_merge_segments(mark.chr, step, overlap, p.thresh, "diff", debug)
        num.mark = num.mark1
	seg.chr =  data.frame(matrix(0, nrow=length(num.mark), ncol=6, dimnames = list(NULL, c("chrom", "loc.start", "loc.end", "num.mark", "length.ratio", "seg.mean"))))
	seg.chr$num.mark = num.mark
	seg.chr$chrom = chr2
	Start = array(1, dim=nrow(seg.chr))
	if (nrow(seg.chr) > 1) {
	    for (i in 2 : nrow(seg.chr)) Start[i] = Start[i - 1] + seg.chr$num.mark[i - 1]
	}
	seg.chr$loc.start = mark.chr$ChrPos[Start] - 49
	seg.chr$loc.end = mark.chr$ChrPos[Start + seg.chr$num.mark - 1] + 50
	for (j in 1:length(Start)) {
	    seg.chr$seg.mean[j] = mean(mark.chr$diff[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	}

	seg.diff = abs(seg.chr$seg.mean[-nrow(seg.chr)] - seg.chr$seg.mean[-1])
	seg.diff.id = which(seg.diff == min(seg.diff))[1]
	while (min(seg.diff) < thresh.diff) {
	    seg.chr$seg.mean[seg.diff.id] = (seg.chr$seg.mean[seg.diff.id] * seg.chr$num.mark[seg.diff.id] + seg.chr$seg.mean[seg.diff.id+1] * seg.chr$num.mark[seg.diff.id+1])/(seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1])
	    seg.chr$num.mark[seg.diff.id] = seg.chr$num.mark[seg.diff.id] + seg.chr$num.mark[seg.diff.id+1]
	    seg.chr$loc.end[seg.diff.id] = seg.chr$loc.end[seg.diff.id + 1]
	    seg.chr = seg.chr[-(seg.diff.id + 1),]
	    seg.diff = seg.diff[-seg.diff.id]
	    if (seg.diff.id > 1)
		seg.diff[seg.diff.id-1] = abs(seg.chr$seg.mean[seg.diff.id] - seg.chr$seg.mean[seg.diff.id-1])
	    if (seg.diff.id <= length(seg.diff))
		seg.diff[seg.diff.id] = abs(seg.chr$seg.mean[seg.diff.id+1] - seg.chr$seg.mean[seg.diff.id])
	    seg.diff.id = which(seg.diff == min(seg.diff))[1]
	}

	seg.chr$length.ratio = seg.chr$num.mark * 100 / (seg.chr$loc.end - seg.chr$loc.start + 1)
	return(seg.chr)
}



########## Main Program ##########

if (exists("norm.chr")) print(norm.chr)
if (exists("norm.seg")) print(norm.seg)
if (exists("est.D")) print(est.D)

######### Set Working Directory ###########

if (exists("working_directory")) {
    setwd(working_directory)
    code.dir = working_directory
} else {
  print("The working directory must be set, Program will abort.")
  q()
}
getwd()

######### Create Location to Store Results ##########

dir.name = "Result"
if (!file.exists(dir.name)) dir.create(dir.name)
paste("Starting time:", date())
quality.merge=F

if (singleBAM) quality.merge=F
if (length(chr.toScale) != length(chr.scale)) {
    print("The length of chr.toScale and chr.scale is different, quit.")
    q()
}

######### Compute Loss of Heterozygosity File ##########

print(paste("SAMPLE", SAMPLE))
file.name = paste(SAMPLE,".ai", sep="")

print(paste("file.name", file.name))

if (!file.exists(file.name)) plotAI=F
if (plotAI) {
    ai = read.table(file.name, header=T)
    ai.file = file.name
    file.name = paste("Result/", SAMPLE, "_LOH_RegTree.txt", sep="")
    if (!file.exists(file.name)) {
	num.test = 0
	offset = array(0, 23)
	
        for (i in 1:22) {
	    num.mark = sum(ai$Chr==paste("chr", i, sep="")) - sum(is.nan(ai$AIDiff[which(ai$Chr==paste("chr", i, sep=""))]))
	    num.test = num.test + num.mark^3/6 - num.mark^2/2 + num.mark/3
	}

	p.thresh = 0.05 / num.test
	step = 20000
	overlap = step/20
	for (i in 1:22) {
	    mark.chr = ai[which(ai$Chr == paste("chr", i, sep="") & (!is.na(ai$AIDiff))),]
	    colnames(mark.chr) = list("Chr", "ChrPos", "diff", "baft", "bafn")
	    current = 1
	    offset = 0
	    offset2 = 0
	    num.mark = -1
	    while (current <= nrow(mark.chr)) {
		if (current + step*1.5 - 1 < nrow(mark.chr)) {
		    mark.tmp = mark.chr[current + ((-offset):(step - 1)),]
	 	    current = current + step
		} else {
		    mark.tmp = mark.chr[(current - offset):length(mark.chr$diff),]
		    current = nrow(mark.chr) + 1
		}
		mark.seg = fun.help(mark.tmp, p.thresh=p.thresh)
		if (num.mark[1] < 0) {
		    num.mark = mark.seg
		} else {
		    mark.seg[1] = mark.seg[1] + offset2
		    if (length(num.mark) > 1) {
			num.mark = c(num.mark[1:(length(num.mark) - 1)], mark.seg)
		    } else {
			num.mark = mark.seg
		    }
		}
		if (num.mark[length(num.mark)] > overlap) {
		    offset = overlap
		} else {
		    offset = num.mark[length(num.mark)]
		}
		offset2 = num.mark[length(num.mark)] - offset
	    }
	    Start = array(1, dim=length(num.mark))
	    if (length(Start) > 1) {
		for (j in 2 : length(num.mark)) Start[j] = Start[j - 1] + num.mark[j - 1]
		p.val = array(1, length(num.mark)-1)
		for (j in 1: length(p.val))
	    	    p.val[j] = t.test(mark.chr$diff[Start[j]:(Start[j] + num.mark[j] - 1)], mark.chr$diff[Start[j + 1]:(Start[j + 1] + num.mark[j + 1] - 1)])$p.value
		p.max.id = which(p.val == max(p.val))[1]
		while (max(p.val) > p.thresh) {
	    	    num.mark[p.max.id] = num.mark[p.max.id] + num.mark[p.max.id+1]
		    num.mark = num.mark[-(p.max.id + 1)]
		    p.val = p.val[-p.max.id]
	    	    Start = Start[-(p.max.id+1)]
	    	    if (p.max.id > 1)
			p.val[p.max.id-1] = t.test(mark.chr$diff[Start[p.max.id-1]:(num.mark[p.max.id - 1] + Start[p.max.id-1] - 1)], mark.chr$diff[Start[p.max.id]:(num.mark[p.max.id] + Start[p.max.id] - 1)])$p.value
	    	    if (p.max.id <= length(p.val))
			p.val[p.max.id] = t.test(mark.chr$diff[Start[p.max.id]:(num.mark[p.max.id] + Start[p.max.id] - 1)], mark.chr$diff[Start[p.max.id+1]:(num.mark[p.max.id+1] + Start[p.max.id+1] - 1)])$p.value
	    	    p.max.id = which(p.val == max(p.val))[1]
		}
	    }
	    seg.chr =  data.frame(matrix(0, nrow=length(num.mark), ncol=5, dimnames = list(NULL, c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean"))))
	    seg.chr$num.mark = num.mark
	    seg.chr$chrom = i
	    Start = array(1, dim=nrow(seg.chr))
	    if (nrow(seg.chr) > 1) {
		for (j in 2 : nrow(seg.chr)) Start[j] = Start[j - 1] + seg.chr$num.mark[j - 1]
	    }
	    seg.chr$loc.start = mark.chr$ChrPos[Start]
	    seg.chr$loc.end = mark.chr$ChrPos[Start + num.mark - 1]
	    for (j in 1 : length(Start)) {
		seg.chr$seg.mean[j] = mean(mark.chr$diff[Start[j]:(Start[j] + seg.chr$num.mark[j] - 1)], na.rm=T)
	    }
	    if (i ==1) {
		seg.final = seg.chr
	    } else {
		seg.final = rbind(seg.final, seg.chr)
	    }
	}
	write.table(seg.final, file.name, col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    }
}

numChr = 22
if (withXY) numChr = 24
offset = array(0, numChr + 1)
norm.all=F

if (singleBAM) {
   norm.all=T
}

if (auto.ref) {
    file.name = paste("Result/", SAMPLE, "_LOH_RegTree.txt", sep="")
    loh = try(read.table(file.name, header=T), silent=T)
    if (class(loh) == "try-error") {
	auto.ref=F
	norm.all = T
    } else {
	id = which(loh$seg.mean < 0.1 & (loh$loc.end - loh$loc.start >= 1e6))
	if (length(id) == 0) {
	    auto.ref=F
	    norm.all = T
	} else {
	    loh = loh[id,]
	    print(loh)
	}
    }
}  else {
    if (exists("norm.chr")) norm.all = (length(norm.chr) == numChr)
}

if (!setGC) {
    if (!exists("gc.prefix")) {
	if (hg18) {
	    gc.prefix = "../../gc/hg18/"
	} else {
    	    gc.prefix = "../../gc/GRCh37/GRCh37-lite-"
	}
	
	if (GRCh38) gc.prefix = "../../gc/GRCh38/GRCh38-"
    }
}

if (!setMap) {
    if (!exists("map.prefix")) {
	if (hg18) {
	    map.prefix = "../../mapability/hg18/"
	} else {
	    map.prefix = "../../mapability/hg19/"
	}
    }
}

softChr = 1
for (chr in 1 : numChr) {
    if (chr < 23) {
     	gc = read.table(paste(gc.prefix, chr, "_", window,".txt", sep=""), header=T)
    }
    if (chr == 23) {
	gc = read.table(paste(gc.prefix, "X_", window,".txt", sep=""), header=T)
    }
    if (chr == 24) {
	gc = read.table(paste(gc.prefix, "Y_", window,".txt", sep=""), header=T)
    }
    file.name<-paste(SAMPLE, "_chr",chr, "_100",sep="")
    if (!file.exists(file.name) & chr==23) file.name<-paste(SAMPLE, "_chrX_100",sep="")
    if (!file.exists(file.name) & chr==24) file.name<-paste(SAMPLE, "_chrY_100",sep="")
    if (singleBAM) {
	file.name<-paste(SAMPLE, "_chr",chr, "_", window,sep="")
 	if (!file.exists(file.name) & chr==23) file.name<-paste(SAMPLE, "_chrX_", window,sep="")
	if (!file.exists(file.name) & chr==24) file.name<-paste(SAMPLE, "_chrY_", window,sep="")
   }
    chr.D =  try(read.table(file.name, header=T), silent=T)
    if (class(chr.D) == "try-error") {
	print(paste("Error in reading", file.name))
	q()
    }
    chr.D = data.frame(Chr = chr, Position = gc$Position + offset[chr], ChrPos = gc$Position, Mean = chr.D[,1], GMean = chr.D[,2], N = gc$N, GC = gc$GC)
    if (mapability) {
	chr.D$Map = scan(paste(map.prefix, "Mapability_chr", chr, "_", window,  sep=""))
    } else {
	chr.D$Map = NA
    }
    chr.D$SoftChr = softChr
    softChr = softChr + 1
    chr.D$Log = log2((chr.D$Mean + 0.01) / (chr.D$GMean + 0.01))
    
    if (sum(chr.D$Mean > 0) == 0) next
    
chr.D = chr.D[which(chr.D$Mean > 0),]
    
    offset[chr + 1] = offset[chr] + nrow(chr.D) * window
    if (chr == 1) {
	cn.D = chr.D
    } else {
	cn.D[nrow(cn.D) + (1:nrow(chr.D)),] = chr.D
    }
    if (auto.ref) {
	id = which(loh$chrom == chr)
	if (length(id) > 0) {
	    for (i in 1:length(id)) {
		id2 = (chr.D$ChrPos >= loh$loc.start[id[i]] & chr.D$ChrPos <= loh$loc.end[id[i]])
		if (length(id2) > 0) {
		    cont = T
		    if (exists("est.D")) {
			if (median(chr.D$Mean[id2], na.rm=T) < est.D * 1.25) cont =F
		    }
		    if (cont) {
			if (!exists("norm.D")) {
			    norm.D = chr.D[id2,]
			} else {
			    norm.D[nrow(norm.D) + (1 : sum(id2)),] = chr.D[id2,]
			}
		    }
		}
	    }
	} 
    } else if (!norm.all) {
	if (exists("norm.seg")) {
	    if (chr == norm.seg[1]) {
		norm.D = chr.D[which(chr.D$ChrPos > norm.seg[2] & chr.D$ChrPos <= norm.seg[3]),]
	    }
	} else {
	    if (chr == norm.chr[1]) {
	    	norm.D = chr.D
	    } else {
	    	if (sum(chr == norm.chr) > 0) {
		    norm.D[nrow(norm.D) + (1:nrow(chr.D)),] = chr.D
	    	}
	    }
	}
    }
    rm(chr.D)
    rm(gc)
    print(paste(chr, softChr - 1, date()))
    if (debug | chr==5 | chr==13) print(gc(debug))
}
cn.orig=cn.D
file.name = ""
paste("Finishing reading genome data:", date())
if (norm.all) {
     m1 = median(cn.D$Mean, na.rm=T)
     cn.D$Mean = cn.D$Mean / m1
     if (!singleBAM) {
	m2 = median(cn.D$GMean, na.rm=T)
	cn.D$GMean = cn.D$GMean / m2
	cn.D$Log = cn.D$Log - median(cn.D$Log, na.rm=T)
     }
} else {
     if (auto.ref) {
	m1 = median(norm.D$Mean, na.rm=T)
       norm.D = norm.D[which(norm.D$Mean <= 1.25 * m1),] 
     }
     m1 = median(norm.D$Mean, na.rm=T)
     cn.D$Mean = cn.D$Mean / median(norm.D$Mean, na.rm=T)
     norm.D$Mean = norm.D$Mean / median(norm.D$Mean, na.rm=T)
     if (!singleBAM) {
     	m2 = median(norm.D$GMean, na.rm=T)
     	cn.D$GMean = cn.D$GMean / median(norm.D$GMean, na.rm=T)
     	cn.D$Log = cn.D$Log - median(norm.D$Log, na.rm=T)
     	norm.D$GMean = norm.D$GMean / median(norm.D$GMean, na.rm=T)
     	norm.D$Log = norm.D$Log - median(norm.D$Log, na.rm=T)
     }
}
if (!singleBAM) {
    cn.D$diff = cn.D$Mean - cn.D$GMean
} else {
    cn.D$diff = cn.D$Mean
}
if (mapability) {
    rm.idx = unique(c(which(cn.D$Map < 0.9), which(is.na(cn.D$GC))))
} else {
  if (gap.norm) {
	rm.tmp = which(cn.D$N > 0.5)
    count = 1
    gap.start = rm.tmp[1]
    current.len = 1
    while (count < length(rm.tmp)) {
	if (rm.tmp[count + 1] - rm.tmp[count] == 1 & (count < length(rm.tmp) - 1)) {
	    if (current.len < gap.thresh) current.len = current.len + 1
	} else {
	    if (gap.start == rm.tmp[1]) {
		rm.idx = max(1, gap.start - current.len) : min(rm.tmp[count] + current.len, nrow(cn.D))
	    } else {
		rm.idx = c(rm.idx, max(1, gap.start - current.len) : min(rm.tmp[count] + current.len, nrow(cn.D)))
	    }
	    gap.start = rm.tmp[count + 1]
	    current.len = 1
	}
	count = count + 1
    }
    rm.idx = unique(rm.idx)
  } else {
    if (singleBAM) {
	rm.idx = unique(c(which(cn.D$N > 0.5), which(cn.D$Mean>4)))
    } else {
	rm.idx = unique(c(which(cn.D$N > 0.5), which(cn.D$Mean>8), which(cn.D$GMean>4)))
    }
  }
}
if (!norm.all) {
  if (mapability) {
    rm.idx2 = unique(c(which(norm.D$Map < 0.9), which(is.na(norm.D$GC))))
  } else {
    if (gap.norm) {
	if (singleBAM) {
	    rm.tmp = which(norm.D$N > 0.5)
	} else {
	    rm.tmp = unique(c(which(norm.D$N > 0.5), which(norm.D$Mean==0 & norm.D$GMean==0)))
	}
    	count = 1
    	gap.start = rm.tmp[1]
    	current.len = 1
    	while (count < length(rm.tmp)) {
	    if (rm.tmp[count + 1] - rm.tmp[count] == 1 & (count < length(rm.tmp) - 1)) {
	    	if (current.len < gap.thresh) current.len = current.len + 1
	    } else {
	    	if (gap.start == rm.tmp[1]) {
		    rm.idx2 = max(1, gap.start - current.len) : min(rm.tmp[count] + current.len, nrow(norm.D))
	    	} else {
		    rm.idx2 = c(rm.idx2, max(1, gap.start - current.len) : min(rm.tmp[count] + current.len, nrow(norm.D)))
	    	}
	    	gap.start = rm.tmp[count + 1]
	    	current.len = 1
	    }
	    count = count + 1
    	}
    } else {
	if (singleBAM) {
	    rm.idx2 = unique(c(which(norm.D$N > 0.5), which(norm.D$Mean>4)))
	} else {
	    rm.idx2 = unique(c(which(norm.D$N > 0.5), which(norm.D$Mean>8), which(norm.D$GMean>4)))
	}
    }
  }
}
cn.D = cn.D[-rm.idx,]
    if (norm.all) {
	cn.lm1 = lm(cn.D$diff~cn.D$GC)
	cn.D$diff = cn.lm1$residuals
	if (!singleBAM) {
	    cn.lm2 = try(lm(cn.D$Log~cn.D$GC), silent=T)
		if (class(cn.lm2) == "try-error") {
		    print("Error in reference selection.  Quitting.")
		    quit(status = 20)
		}
	 	cn.D$Log = cn.lm2$residuals
	}
    } else {
	gc = norm.D$GC
	if (singleBAM) {
	    cn.lm1 = try(lm((norm.D$Mean)~gc, subset=-rm.idx2), silent=T)
	} else {
	    cn.lm1 = try(lm((norm.D$Mean-norm.D$GMean)~gc, subset=-rm.idx2))
	}
        if (class(cn.lm1) == "try-error") {
            print("Error in reference selection.  Quitting.")
            quit(status = 20)
        }

	pre.tmp = data.frame(gc = cn.D$GC)
 	cn.D$diff = cn.D$diff - predict(cn.lm1, pre.tmp)
	if (!singleBAM) {
 		cn.lm2 = try(lm((norm.D$Log)~gc, subset=-rm.idx2), silent=T)
            if (class(cn.lm2) == "try-error") {
                print("Error in reference selection.  Quitting.")
                quit(status = 20)
            }

 	    cn.D$Log = cn.D$Log - predict(cn.lm2, pre.tmp)
	}
	rm(pre.tmp)
    }
rm(cn.lm1)
if (!singleBAM) rm(cn.lm2)
gc(T)
paste("Finishing correcting GC content for difference:", date())
if (!setPThresh) {
    num.test = 0
    for (i in 1 : numChr) {
	num.mark = sum(cn.D$Chr==i)
	num.test = num.test + num.mark^3/6 - num.mark^2/2 + num.mark/3
    }
    p.thresh = 0.05/num.test
}
if (singleBAM) {
    file.prefix = paste(dir.name,"/",SAMPLE, "_Single","_CONSERTING_", sep="")
} else {
    file.prefix = paste(dir.name,"/",SAMPLE, "_CONSERTING_", sep="")
}
if (debug) file.prefix = paste(file.prefix, "Debug_", sep="") 
if (mapability) {
	file.prefix = paste(file.prefix, "Mapability_",window, sep="")
} else {
	file.prefix = paste(file.prefix,window, sep="")
}    
file.name = paste(file.prefix, ".txt",sep="")

if (!merge.diff) file.name = paste(file.prefix, "_raw.txt",sep="")
file.final=""
if (forced | (!file.exists(file.name))) {
    seg.final=NULL
    for (i in 1:numChr) {
       print(paste("Chr",i,":    ",date()))
       if (debug) print(gc(T))
       if (singleBAM)  {
           seg.chr = fun.chr.germline(i, p.thresh=p.thresh)
       } else {
           seg.chr = fun.chr(i, p.thresh=p.thresh)
       }
       if (i ==1) {
          seg.final = seg.chr
       } else {
          seg.final = rbind(seg.final, seg.chr)
       }
    }

    seg.final[,-(1:4)] = round(seg.final[,-(1:4)] * 1000) / 1000
    
     drops <- c("length.ratio")
     tmp.final <- seg.final[ , !(names(seg.final) %in% drops)]
     seg.final <- tmp.final

    write.table(seg.final, file.name, col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    seg.final$chrom[which(seg.final$chrom==23)] = "X"
    seg.final$chrom[which(seg.final$chrom==24)] = "Y"
}
gc(T)

seg.final = read.table(file.name, header=T)
if (!singleBAM) {
    id.Y = which(seg.final$chrom==24)
    if (sum(seg.final$num.mark[id.Y] * seg.final$GMean[id.Y])/sum(seg.final$num.mark[id.Y]) > 0.2) Male = T
}

if (quality.merge) {
    if (Male) {
	system(paste("java.sh -cp ", code.dir, " CNVQualityMerge -c1", sv.file, "-seg", file.name, "-c2", file.final, "-ai", ai.file, "-male -checkDmeans", "-td", thresh.diff, "-tl", logthresh.diff, "-ts", ts))
    } else {
	system(paste("java.sh -cp ", code.dir, " CNVQualityMerge -c1", sv.file, "-seg", file.name, "-c2", file.final, "-ai", ai.file, " -checkDmeans", "-td", thresh.diff, "-tl", logthresh.diff, "-tl", logthresh.diff, "-ts", ts))
    }
    file.name = paste(file.name, ts, ".QualityMerge", sep="")
    system(paste("java.sh -cp ", code.dir, " ReformatGeDI ", file.name))
}
seg.final = read.table(file.name, header=T)
paste("Finishing building tree:", date())
print(paste("Job completed:", date()))
proc.time()
