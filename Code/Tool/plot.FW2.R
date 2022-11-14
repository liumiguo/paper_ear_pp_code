#====================================================================================
#plot.FW2
#====================================================================================

#按照FW包中的plot.FW函数进行修改的函数
#要求plotVAR中将野生型放在最后面

#x 为FW模型拟合的对象
#plotVAR 要展示的VAR名称，但是野生型编号要放在最后面

# EarLength2018_Genotype<-EarLength_g2018%>%
#   filter(ConstructID%in%c(Construct_ID,"99999"), G>=.[ConstructID=="99999",]$G)
# plot.FW2(x=EarLength2018,plotVAR=c(EarLength2018_Genotype$Genotype),main=str_c(Construct_ID,"\n"))


plot.FW2<-function(x,chain=1,plotVAR,ENVlabel = FALSE,...){
  # x<-EarLength2018
  # chain = 1
  # ENVlabel="split"
  # plotVAR<-c("0233321009","9111421013","9999901001")
  y = x$y#y实测值
  VARlevels = x$VARlevels#有多少个var就多少个标签
  ENVlevels = x$ENVlevels#有多少个ENV就有多少个标签
  VAR = factor(x$VAR, levels = VARlevels)#创建VAR并使用x$VARlevels排序
  ENV = factor(x$ENV, levels = ENVlevels)
  IDL = as.numeric(VAR)#转化为数字化
  IDE = as.numeric(ENV)
  yhat = x$yhat[, chain]#y预测值
  mu = x$mu[chain]
  g = x$g[, chain]#g值
  names(g) = VARlevels
  b = x$b[, chain]
  names(b) = VARlevels
  h = x$h[, chain]
  names(h) = ENVlevels



  if (!is.null(plotVAR)) {#将字符串转化为数字
    if (is.integer(plotVAR)) {
      plotIDL = plotVAR
    }
    else if (is.character(plotVAR)) {
      plotIDL = match(plotVAR, VARlevels, nomatch = 0)
    }
  }else {#如果空则如此
    plotIDL = c(1:length(VARlevels))
  }

  n.plotVAR = length(plotIDL)
  oripar = par()$mar
  par(mar = oripar)
  if (ENVlabel == FALSE) {#为是否显示ENVlabel存放
    par(xpd = T, mar = oripar + c(1.5, 0, 0, 5))
  }else {
    par(xpd = T, mar = oripar + c(0, 0, 0, 5))
  }

  range.h = range(h, na.rm = T)#找环境数值的区间
  y.plot = y[(IDL %in% plotIDL)]#确定哪些品种要绘制的y值
  range.y = range(y.plot, na.rm = T)#确定y值的区间
  args1 = list(xlab = "Environment effect", ylab = "Variety performance",
               type = "n")#指定绘图的标签
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  # do.call(plot, c(list(formula = range.y ~ range.h), args1))
  if (ENVlabel == FALSE) {#设置图形的x,y轴信息
    do.call(plot, c(list(formula = range.y ~ range.h), args1))
  }else {
    args2 = args1
    args2$xlab = ""
    do.call(plot, c(list(formula = range.y ~ range.h), args2))
    args2 = args1
    args2$main = ""
    args2$ylab = ""
    args2$line = 4
    args2$type = NULL
    do.call(title, args2)
  }

  sorth = sort(h)#环境从小到大排列
  sorth1 = sorth[seq(1, length(h), by = 2)]#奇数环境
  sorth2 = sorth[seq(2, length(h), by = 2)]#偶数环境

  # if (ENVlabel != FALSE) {#添加环境的上下标签
  #   if (ENVlabel == "split" | ENVlabel == TRUE) {
  #     axis(side = 1, at = sorth1, labels = names(sorth1),
  #          line = 2)
  #     axis(side = 3, at = sorth2, labels = names(sorth2),
  #          line = 1)
  #   }
  #   else if (ENVlabel == "bottom") {
  #     axis(side = 1, at = sorth, labels = names(sorth),
  #          line = 2)
  #   }
  #   else if (ENVlabel == "top") {
  #     axis(side = 3, at = sorth, labels = names(sorth),
  #          line = 1)
  #   }
  #   else {
  #     cat("ENVlabel must be TRUE, FALSE, \"split\",\"top\" or \"bottom\" \n")
  #   }
  # }

  cols = NULL
  pchs = NULL
  usr <- par("usr")
  col = 1

  #------
  for (i in plotIDL) {#plotIDL是所绘制材料在第几个
    # i<-plotIDL
    col = col + 1
    pch = 1
    cols = c(cols, col)
    pchs = c(pchs, pch)
    IDLi = which(IDL == i)#IDL为将VAR因子数字化的，该行为将选择指定的VAR所在的行
    y.i = y[IDLi]#筛选指定的VAR所在的行的y值（实测值）
    ENV.i = ENV[IDLi]#筛选指定的VAR所在的行的ENV值（实测值）
    cellMeans = aggregate(y.i, by = list(ENV.i), mean)
    cellMeansIDE = as.numeric(cellMeans[, 1])
    cellMeans = cellMeans[, 2]
    hi = h[cellMeansIDE]
    clip(range.h[1], range.h[2], range.y[1], range.y[2])
    if(i==plotIDL[length(plotIDL)]){
      abline(a = mu + g[i], b = b[i] + 1, col = 1,lwd = 2,...)
    }else{
      abline(a = mu + g[i], b = b[i] + 1,lwd = 2, col = col,...)#!!!!!!!!!!!!!!!!!
    }
    do.call("clip", as.list(usr))
    if(i==plotIDL[length(plotIDL)]){
      points(cellMeans ~ hi, col = "black", pch = pch, xpd = T,...)#!!!!!!!!!!!!!!
    }else{
      points(cellMeans ~ hi, col = col, pch = pch, xpd = T,...)#!!!!!!!!!!!!!!
    }

  }
  clip(range.h[1], range.h[2], range.y[1], range.y[2])
  # abline(a = mean(y.plot, na.rm = T), b = 1, lty = 2, col = 1)
  clip(usr[1], usr[2] + (usr[2] - usr[1]) * 10, usr[3], usr[4])
  legend(x = usr[2], y = usr[4], legend = c(VARlevels[plotIDL]),
         lty = c(rep(1, n.plotVAR)),
         col = c(cols[1:length(cols)-1],1), bg = "transparent", bty = "n")
  par(mar = oripar)
  do.call("clip", as.list(usr))
}
