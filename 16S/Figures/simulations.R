#!/usr/bin/env Rscript

#suppressMessages(library("optparse"))
#option_list = list(make_option(c("-v","--value"), action="store", type='character', 
#                              help="simulated value, corresponding to 'value_dir', part of directory structure [required]"),
#                  make_option(c("-s","--seed"), action="store", type='character', 
#                              help="seed(s, dot-sep) used in simulation, corresponding to 'seed_dir', part of directory structure: 'value_dir/sim_dir' [required]"),
#                  make_option(c("-V","--vddir"), action="store", type='character', 
#                              help="path to directory with VD results, parent dir of 'value_dir' [required]"),
#                  make_option(c("-D","--rdata"), action="store", type='character', 
#                              help="name of rdata in which saved estimates and STE, [required]"),
#                  make_option(c("-S","--simdir"), action="store", type='character', 
#                              help="path to directory with simulations, parent dir of 'value_dir' [required]"),
#                  make_option(c("-O","--out"), action="store", type='character', 
#                              help="path to file with objects for plotting [required]"),
#                  make_option(c("-m","--model"), action="store", default = "bi", type='character', 
#                              help="model used in VD, 'bi' or 'uni'[TODO check if working], with '_null' if null [default: %default]") # TODO: check it works with univariate too
#)

sourcefun = "/users/abaud/htonnele/git/me/core_VD/code/Rsrc/functions/" # 
source(file.path(sourcefun, "prepare_res.R"))

##### MI - different corr ######
opt=list(pvar="cor(DGE,IGE)",
         value="0.9,0.0,neg0.9", 
         pop = "MI",
         seed= "21",
         vddir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_422/setUhDGhIG/DG1_IG1/",
         simdir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_422/setUhDGhIG/DG1_IG1/",
         outpre = "/users/abaud/htonnele/PRJs/plots/simP50/setUhDGhIG/paper/MI_DG1_IG1",
         model = "uni")

##### NY - different corr ######
#opt=list(pvar="cor(DGE,IGE)",
#         value="0.9,0.0,neg0.9",
#         pop = "NY",
#         seed= "22",
#         vddir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_654/setUhDGhIG/DG1_IG1/",
#         simdir = "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_654/setUhDGhIG/DG1_IG1/",
#         outpre = "/users/abaud/htonnele/PRJs/plots/simP50/setUhDGhIG/paper/NY_DG1_IG1",
#         model = "uni")#"uni")

####### get options for all #####
val = unlist(strsplit(opt$value,","))
pdfVCs = paste0(opt$outpre, "_VCs_from_sim_", paste(val,collapse="."),".pdf")
#pdfTotV = paste0(opt$outpre, "_totV_from_sim_",paste(val,collapse="."),".pdf")

pvar=opt$pvar
pop=opt$pop
sid = unlist(strsplit(opt$seed,"[.]"))   # Otherwise, keep as character
vddir = opt$vddir
simdir = opt$simdir
model = unlist(strsplit(opt$model, "_"))[1]


###### VD files - different IGE vs DGE ######
get_res = function(vddir, val, sid, simdir){
  nocols=c(grep("2",colest,value=T),"covariates_names","conv","LML","time_exec","union_focal", "inter_focal", "union_cm", "inter_cm")
  
  get_est = function(namerdata){
    estE = data.frame()
    for(v in val){
      ## getting VD results
      f = list.files(file.path(vddir, v, sid), pattern = namerdata, 
                           recursive=T,full.names=T)
      load(f);
      if(exists("res")){
        VCs=res$VCs
      }
      est_f = VCs[,! colnames(VCs) %in% nocols]
      est_f[,"simcor"] = rep(gsub("neg","-",v), nrow(est_f))
      # binding to the big dataframe
      estE = rbind(estE, est_f)
      
      rm(est_f, VCs)
    }
    return(estE)
  }
  
  DGEest = get_est("P50_Rn7_pruned_DGE_cageEffect_None_all_estNste.Rdata")
  IGEest = get_est("P50_Rn7_pruned_DGE_IGE_cageEffect_None_all_estNste.Rdata")
  
  x = IGEest
  y = DGEest
  all_res = data.frame("analysis"=c(rep("wtIGE", nrow(x)), rep("DGEonly", nrow(y))), 
                       "DGE" = c(x[,"prop_Ad1"], y[,"prop_Ad1"]),
                       "IGE" = c(x[,"prop_As1"], rep(NA, nrow(y))), 
                       "cor.DGE.IGE" = c(x[,"corr_Ad1s1"], rep(NA, nrow(y))),
                       
                       "DEE" = c(x[,"prop_Ed1"], y[,"prop_Ed1"]),
                       "IEE" = c(x[,"prop_Es1"], rep(NA, nrow(y))), 
                       "cor.DEE.IEE" = c(x[,"corr_Ed1s1"], rep(NA, nrow(y))),
                       
                       "CE" = c(x[,"prop_C1"], y[,"prop_C1"]), 
                       "tot.phenot.var" = c(x[,"total_var1"], y[,"total_var1"]),
                       "simcor" = c(x[,"simcor"], y[,"simcor"]))
  
    
    all_res[,"analysis"] = factor(all_res[,"analysis"], levels=c("DGEonly", "wtIGE"))
    all_res[,"simcor"] = factor(all_res[,"simcor"], levels=unique(all_res[,"simcor"])[order(unique(all_res[,"simcor"]))])#c("-0.9","0.0","0.9"))
    
    ##### Simulated values #### 
    prop_names = grep("prop_", colnames(IGEest), value = T)
    corr_names = grep("corr_", colnames(IGEest), value = T) # NB: this is equal to corr_uni_names when have estimates from uni
    totv_names = grep("tot_|total_", colnames(IGEest), value = T)
    corP_names = grep("corParams", colnames(IGEest), value=T)
    
    sim_param = data.frame()
    for(v in val){ # need IGEest - or estimate of some sort 
      simfile = list.files(file.path(simdir, v, sid), pattern = "params", recursive =T, full.names=T)
      
      #b. get simulated params -> put in sim_prop, sim_corr, sim_tot
      simP_df = read.csv(simfile, header = T, sep="\t")
      sim_prop = simP_df[gsub("prop_", "var_", prop_names), "prop_params"]
      names(sim_prop) = prop_names
      
      sim_corr = simP_df[corr_names, "set_params"] # not splitting yet in corr_uni, corr_bi as dataframe from uni has only phenotype 1
      names(sim_corr) = corr_names 
      
      sim_totv = simP_df[c("var_y1"), "var_term"] # NB: this can't compare as the rest is stdized! ??
      names(sim_totv) = c("total_var1") #,"total_var2" ) #totv_names
      
      sim_param = rbind(sim_param,
                        data.frame("param" = names(c(sim_prop, sim_corr, sim_totv)),
                                   "values" = c(sim_prop, sim_corr, sim_totv), 
                                   "corr" = rep(gsub("neg","-",v), length(c(sim_prop, sim_corr, sim_totv)))))
    }
    #str(sim_param)
    return(list("all_res" = all_res, "sim_param" = sim_param))
}

    
# defining ylim for each plot - so that concordant between the two populations
all_resMI = get_res("/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_422/setUhDGhIG/DG1_IG1/", val, 21, 
                    "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_422/setUhDGhIG/DG1_IG1/")
all_resNY = get_res("/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_654/setUhDGhIG/DG1_IG1/", val, 22,
                    "/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_654/setUhDGhIG/DG1_IG1/")
ylimi = sapply(c("DGE","IGE","cor.DGE.IGE","DEE","IEE","cor.DEE.IEE","CE","tot.phenot.var"), 
               \(p) range(c(all_resMI$all_res[,p], all_resNY$all_res[,p]), na.rm = T))
rm(all_resMI, all_resNY)
#ylimi = data.frame("DGE"= c(0.0672384,0.6089855),        
#                   "IGE"= c(0.05572387,0.36439842), 
#                   "cor.DGE.IGE"= c(-0.9996709,0.9996614),
#                   "DEE"=c(9.717279e-05, 7.401876e-01 ),          
#                   "IEE"=c(9.060918e-05, 4.420130e-01), 
#                   "cor.DEE.IEE"= c(-0.9988256, 0.9984923),           
#                   "CE"= c(4.153922e-09,7.550063e-01),      
#                   "totv"=c(0.8576124,1.1543911))

res = get_res(vddir, val, sid, simdir)
all_res = res$all_res
ylimi= sapply(c("DGE","IGE","cor.DGE.IGE","DEE","IEE","cor.DEE.IEE","CE","tot.phenot.var"), \(p) range(all_res[,p], na.rm = T))
sim_param = res$sim_param
head(all_res) # check
table(all_res[,c("analysis","simcor")]) # check
nsim = unique(table(all_res[,c("analysis","simcor")]))



#### ... to plot #####
#boxcol = c("analysed with DGE only"="#4DBBD5FF", "analysed with DGE and IGE"="#FFAA0EFF")
boxcol = c("analysed with DGE only"="#3C5488FF", "analysed with DGE and IGE"="#00A087FF")
#segcol = c("#DC0000FF")
segcol=c("#DC0000FF")
ats = data.frame("at" = c(1:2,4:5,7:8),
                 "analysis" = rep(c("DGEonly", "wtIGE"), 3),
                 #"simcor" = c(rep("0.01",2), rep("0.1",2), rep("0.4",2)) )  
                 "simcor" = c(sapply(levels(all_res$simcor), \(l) rep(l, 2))))  

#### Function to plot when wtIGE and DGE only ######
complot= function(p,psim,pvar, lg.pos=NULL){
  formula = formula(all_res[,p] ~ all_res[,"analysis"] + all_res[,"simcor"])
  #ats = c(1:2, 4:5, 7:8) #where to put boxplots
  ylim = ylimi[,p]
  vioplot::vioplot(formula, col = alpha(boxcol,0.2), #col = rep(pal, 3),
                   at = ats[,"at"], names = rep('', length(ats[,"at"])), 
                   las = 1, xlab = '', xaxt = "n", ylab='',
                   cex.axis = 1.25,  
                   border = "black", 
                   drawRect=F,
                   wex=0.5,
                   areaEqual = F,
                   ylim = ylim)
  abline(h=mean(sim_param[sim_param$param==psim,"values"]), col="grey", lty=1, lwd=0.5)
  xcoord = c(sapply(ats[,"at"], function(a){rep(a, nsim)}))
  points(x=jitter(xcoord, factor=0.2), 
         y=c(all_res[order(all_res$simcor,all_res$analysis),p]), cex=0.5, col=rep(sapply(boxcol, \(a) {rep(a,nsim)}), 3) )
  boxplot(formula, at=ats[,"at"], xlim=c(0,9), 
          add=T, 
          boxwex = 0.2, col=alpha("white",0.5), outline=F, 
          xaxt="n", yaxt="n")
  segments(x0=ats[c(1,3,5),"at"]-0.5, 
           y0=sim_param[order(sim_param$corr),][sim_param$param==psim,"values"],
           x1=ats[c(2,4,6),"at"]+0.5,
           y1=sim_param[order(sim_param$corr),][sim_param$param==psim,"values"], 
           col=segcol,
           #col=c("#E64B35FF", "#FFAA0EFF", "#FFD200FF"), 
           lty=2, lwd=2)
  axis(side = 1, at = c(sapply(unique(ats[,"simcor"]), \(c){mean(ats[ats$simcor==c,"at"])})), 
       labels = unique(ats[,"simcor"]), tick = T, cex.axis=1.25)
  title(ylab = paste0("estimated ", p),       
        cex.lab = 1.4,
        line = 3.5)
  title(xlab=paste0("simulated ",pvar),
        cex.lab = 1.4)
        #line = 2.5)
  title(sub = pop) #paste0(nsim, " simulations - ", pop))
  if(is.null(lg.pos)){
    xc = grconvertX(0.22, "npc")  # 0.22 is because it looks aligned on the right with pdf(w=7,h=6) # Use normalized plot coordinates
    yc = grconvertY(1, "nic")    # Position at 10% from bottom
    
    legend(x=xc, y=yc,
           legend=c(paste0("simulated ",p)),
           #fill=c(boxcol,NA), 
           lty=c(2), lwd=c(2), col = c(segcol), seg.len = 1.5,
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           bty="n", x.intersp = 0.5)
    legend(x=(5.2/xc)*xc, y=yc, # 5.2 is because it looks aligned on the right with pdf(w=7,h=6)
           legend=c(names(boxcol)),
           fill=c(boxcol), 
           #lty=c(NA,NA,2), lwd=c(NA,NA,2), col = c(NA,NA,segcol),
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           bty="n", x.intersp = 0.8, horiz = T)
  }else{
    legend(lg.pos, #adj=0, #xjust=0,
           legend=c(paste0("simulated ",p), names(boxcol)),
           #fill=c(),
           text.col = c("black", rep(alpha("white", 0),2)),
           lty=c(2,NA,NA), lwd=c(2,NA,NA), col = c(segcol,NA,NA), seg.len = 1.5,
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           #bty="n", 
           x.intersp = 0.5)
           
    legend(lg.pos, #adj=0,
           legend=c(NA,names(boxcol)),
           fill=c(NA, boxcol), 
           cex=1, border = F, ncol=1, xpd=T, #bg="white",
           bty="n", x.intersp = 1.2)
  }
  
  
  ##legend("topleft",
  ##       legend=c(names(boxcol), paste0("simulated ",p)),
  ##       fill=c(boxcol,NA), 
  ##       lty=c(NA,NA,2), lwd=c(NA,NA,2), col = c(NA,NA,segcol),
  ##       border = F, ncol=2, xpd=T, bg="white", 
  ##       inset = c(0,-0.09))
  #legend(x=xc-0.4*xc,y=yc-0.06*yc,
  #       legend=paste0("simulated ",p),
  #       lty=2, lwd=2, col = segcol,
  #       ncol=1, bty="n", xpd=T)
}

#### Function to plot when only wtIGE  ######
complot2 = function(p,psim,pvar, lg.pos=NULL){
  #ats = c(1:2, 4:5, 7:8) #where to put boxplots
  ylim = ylimi[,p]
  ats = ats[ats$analysis =="wtIGE",]
  boxcol = boxcol["analysed with DGE and IGE"]
  all_res = all_res[all_res$analysis =="wtIGE",]
  
  formula = formula(all_res[,p] ~ all_res[,"simcor"])
  vioplot::vioplot(formula, col = alpha(boxcol,0.3), #col = rep(pal, 3),
                   at = ats[,"at"], 
                   las = 1, xlab = '', xaxt = "n", ylab='',
                   cex.axis = 1.25,  
                   border = "black", 
                   drawRect=F,
                   wex=0.5,
                   areaEqual = F,
                   ylim = ylim)
  #abline(h=mean(sim_param[sim_param$param==psim,"values"]), col="grey", lty=1, lwd=0.5)
  xcoord = c(sapply(ats[,"at"], function(a){rep(a, nsim)}))
  points(x=jitter(xcoord, factor=0.2), 
         y=c(all_res[order(all_res$simcor,all_res$analysis),p]), cex=0.5, col=rep(sapply(boxcol, \(a) {rep(a,nsim)}), 3) )
  boxplot(formula, at=ats[,"at"], xlim=c(0,9), 
          add=T, 
          boxwex = 0.2, col=alpha("white",0.5), outline=F, 
          xaxt="n", yaxt="n")
  segments(x0=ats[,"at"]-0.5, 
           y0=sim_param[order(sim_param$corr),][sim_param$param==psim,"values"],
           x1=ats[,"at"]+0.5,
           y1=sim_param[order(sim_param$corr),][sim_param$param==psim,"values"], 
           col=segcol,
           #col=c("#E64B35FF", "#FFAA0EFF", "#FFD200FF"), 
           lty=2, lwd=2)
  axis(side = 1, at = c(sapply(unique(ats[,"simcor"]), \(c){mean(ats[ats$simcor==c,"at"])})), 
       labels = unique(ats[,"simcor"]), tick = T, cex.axis=1.25)
  title(ylab = paste0("estimated ", p),       cex.lab = 1.4,
        line = 3.5)
  title(xlab=paste0("simulated ",pvar),
        cex.lab = 1.4)#,
        #line = 2.5)
  title(sub = pop) #paste0(nsim, " simulations - ", pop))
  if(is.null(lg.pos)){
    xc = grconvertX(0.22, "npc")  # 0.22 is because it looks aligned on the right with pdf(w=7,h=6) # Use normalized plot coordinates
    yc = grconvertY(1, "nic")    # Position at 10% from bottom
    
    legend(x=xc, y=yc,
           legend=c(paste0("simulated ",p)),
           #fill=c(boxcol,NA), 
           lty=c(2), lwd=c(2), col = c(segcol), seg.len = 1.5,
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           bty="n", x.intersp = 0.5)
    legend(x=(5.2/xc)*xc, y=yc, # 5.2 is because it looks aligned on the right with pdf(w=7,h=6)
           legend=c(names(boxcol)),
           fill=c(boxcol), 
           #lty=c(NA,NA,2), lwd=c(NA,NA,2), col = c(NA,NA,segcol),
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           bty="n", x.intersp = 0.8, horiz = T)
  }else{
    legend(lg.pos, #adj=0, #xjust=0,
           legend=c(paste0("simulated ",p), names(boxcol)),
           #fill=c(),
           text.col = c("black", rep(alpha("white", 0),2)),
           lty=c(2,NA,NA), lwd=c(2,NA,NA), col = c(segcol,NA,NA), seg.len = 1.5,
           cex=1, border = F, ncol=1, xpd=T, bg="white",
           #bty="n", 
           x.intersp = 0.5)
    
    legend(lg.pos, #adj=0,
           legend=c(NA,names(boxcol)),
           fill=c(NA, boxcol), 
           cex=1, border = F, ncol=1, xpd=T, #bg="white",
           bty="n", x.intersp = 1.2)
  }
  
  ##xc = grconvertX(0, "npc")  # Use normalized plot coordinates
  ##yc = grconvertY(0.95, "nic")    # Position at 10% from bottom
  ##legend(x=xc, y=yc,
  ##       legend=c(names(boxcol)),
  ##       fill=c(boxcol), 
  ##       cex = 1.2,
  ##       #lty=c(NA,NA,2), lwd=c(NA,NA,2), col = c(NA,NA,segcol),
  ##       border = F, ncol=1, xpd=T, bg="white",
  ##       bty="n", x.intersp = 1)
  ###inset = c(0,-0.09))
  ##legend(x=7*xc, y=yc,
  ##       legend=c(paste0("simulated ",p)),
  ##       #fill=c(boxcol,NA), 
  ##       lty=c(2), lwd=c(2), col = c(segcol),
  ##       border = F, ncol=1, xpd=T, bg="white",
  ##       cex = 1.2,
  ##       bty="n",
  ##       x.intersp = 0.5)
  ##legend("topleft",
  ##       legend=c(names(boxcol), paste0("simulated ",p)),
  ##       fill=c(boxcol,NA), 
  ##       lty=c(NA,NA,2), lwd=c(NA,NA,2), col = c(NA,NA,segcol),
  ##       border = F, ncol=2, xpd=T, bg="white", 
  ##       inset = c(0,-0.09))
  #legend(x=xc-0.4*xc,y=yc-0.06*yc,
  #       legend=paste0("simulated ",p),
  #       lty=2, lwd=2, col = segcol,
  #       ncol=1, bty="n", xpd=T)
}


pdf(pdfVCs, h=6, w=7)

par(mar=c(5.1,5.1,2.1,2.1))
complot("DGE", "prop_Ad1",pvar, "bottomleft")
complot("DEE", "prop_Ed1",pvar, "topright")
complot("CE", "prop_C1",pvar, "topleft")

###### no DGE only
complot2("IGE", "prop_As1", pvar, "bottomleft") # "bottomleft" #with NY # "topleft" with MI
complot2("IEE", "prop_Es1", pvar, "topleft")
complot2("cor.DGE.IGE", "corr_Ad1s1", pvar, "topleft")
complot2("cor.DEE.IEE", "corr_Ed1s1", pvar, "topleft")
#dev.off()

#pdf(pdfTotV, w=12,h=8)
#par(mar=c(5.1,5.1,4.1,2.1))
complot("tot.phenot.var","total_var1", pvar, "topright")
dev.off()

#the legend is sometimes overlapping the dots (see tot var and cor.DEE.IEE in NY and cor.DEE.IEE in MY) 
#- not the best -  do you think still acceptable?

######### Total var #########
#### y = res_DGE[,"total_var1"]; x= res_IGE[,"total_var1"]
#### xlab="total phenot var (with IGE)"
#### ylab="total phenot var (DGE alone)"
#### xylimi = range(c(x,y))
#### plot(x,y,
####      ylim=xylimi, xlim=xylimi, 
####      xlab=xlab, ylab=ylab)
#### abline(a=0, b=1, col="grey", lty=2)
#### abline(h=sim_totv, v=sim_totv, lty=2, lwd=2, col="red")
#### legend(x = 'topleft',  legend = "simulated total phenot var", lty=2, lwd=2, col="red", cex = 1, bg="white", bty="n")
#### cor.test(x,y)
