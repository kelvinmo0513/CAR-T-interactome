#!/usr/bin/env Rscript
# Render newly computed manuscript tables using the approved compact visual style.
# Rendering never fits a model or substitutes archived numerical results.
suppressPackageStartupMessages({library(grid); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript render.R NUMERICAL_DIR OUTPUT_DIR")
inp <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = FALSE)
if (inp == out) stop("Numerical and plot output directories must differ")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Plot output directory must be empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
qa <- file.path(out, "audit")
dir.create(qa)
read <- function(n) read.csv(file.path(inp, n), stringsAsFactors = FALSE)
stats <- read("model_statistics.csv")
pools <- read("pooled_statistics.csv")
weights <- read("pool_weights.csv")
km <- read("km_curves.csv")
risk <- read("risk_tables.csv")
ks <- read("KM_statistics.csv")
stopifnot(nrow(ks) == 4, all(ks$cohort == "li_cart"), all(ks$n == 22),
          all(stats$estimable), !anyDuplicated(stats$cohort),
          !anyDuplicated(pools$pool))
labels<-c(CCL8='CCL8+',CCL13='CCL13+',CCL8_CCL13='CCL8+CCL13+',CCL18='CCL18+')
colors<-c(High='#FC0001',Low='#0000FC')
txt<-function(s,x,y,size=6.3,bold=FALSE,just='left',col='black',rot=0)
 grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,rot=rot,
 gp=gpar(fontfamily='Helvetica',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
ln<-function(x,y,col='black',lwd=.5,lty=1)
 grid.lines(unit(x,'mm'),unit(y,'mm'),gp=gpar(col=col,lwd=lwd,lty=lty))
box<-function(L,R,B,T)grid.rect(unit(L,'mm'),unit(B,'mm'),unit(R-L,'mm'),unit(T-B,'mm'),
 just=c('left','bottom'),gp=gpar(fill=NA,col='black',lwd=.65))
measure<-function(id){a<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE)
 b<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
 list(id=id,bbox_pt=c(a$x,a$y,b$x,b$y)*72)}
layout_record<-function(name,w,h,panels,columns=list(),rows=list()) {
  # Measure the final device coordinates, not pre-layout guesses. This audit is
  # portable and independent of any installed authoring helper package.
  manifest <- list(schema_version=1,backend='R native grid',
    figure=list(width_pt=w/25.4*72,height_pt=h/25.4*72),
    panels=panels,column_groups=columns,row_groups=rows)
  get_box <- function(id) panels[[match(id, vapply(panels, function(p) p$id, character(1)))]]$bbox_pt
  residuals <- numeric()
  for (g in columns) {
    boxes <- do.call(rbind, lapply(g$panels, get_box))
    residuals <- c(residuals, diff(range(boxes[,1])), diff(range(boxes[,3])))
  }
  for (g in rows) {
    boxes <- do.call(rbind, lapply(g$panels, get_box))
    residuals <- c(residuals, diff(range(boxes[,2])), diff(range(boxes[,4])),
                   diff(range(boxes[,3]-boxes[,1])))
  }
  deviation <- if (length(residuals)) max(residuals) else 0
  if (deviation > 1.5) stop("Rendered panel alignment exceeds 1.5 pt")
  write_json(manifest,file.path(qa,paste0(name,'.alignment-layout.json')),
             auto_unbox=TRUE,pretty=TRUE,digits=12)
  write_json(list(status='PASS',maximum_deviation_pt=deviation,tolerance_pt=1.5),
    file.path(qa,paste0(name,'.alignment.json')),auto_unbox=TRUE,pretty=TRUE,digits=12)
}
export<-function(name,w,h,draw){
 # Base PDF is required because the local Cairo device cannot load libSM.
 pdf(file.path(out,paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Helvetica',useDingbats=FALSE)
 grid.newpage();draw(TRUE);dev.off()
 svglite::svglite(file.path(out,paste0(name,'.svg')),width=w/25.4,height=h/25.4)
 grid.newpage();draw(FALSE);dev.off()
 ragg::agg_png(file.path(out,paste0(name,'.png')),width=w,height=h,units='mm',res=600,background='white')
 grid.newpage();draw(FALSE);dev.off()
 ragg::agg_png(file.path(qa,paste0(name,'_preview.png')),width=w,height=h,units='mm',res=180,background='white')
 grid.newpage();draw(FALSE);dev.off()
}
fmtp<-function(p)format(signif(p,3),trim=TRUE,scientific=p<.001)
draw_km<-function(marker,id=marker){
 s<-ks[ks$cohort=='li_cart' & ks$marker==marker,]
 z<-km[km$cohort=='li_cart' & km$marker==marker,]
 r<-risk[risk$cohort=='li_cart' & risk$marker==marker,]
 txt(labels[[marker]],30,68.5,8,TRUE,'centre')
 for(g in c('High','Low')){
   x<-if(g=='High')15.2 else 34
   ln(c(x,x+3),c(64.5,64.5),colors[[g]],1)
   ln(rep(x+1.5,2),64.5+c(-.5,.5),colors[[g]],.6)
   txt(g,x+4.5,64.5,6.5,col=colors[[g]])
 }
 L<-11;R<-58;B<-28.5;T<-60
 xr<-c(-.045*42,1.02*42);yr<-c(-.04,1.04)
 xp<-function(t)L+(t-xr[1])/diff(xr)*(R-L)
 yp<-function(v)B+(v-yr[1])/diff(yr)*(T-B)
 ticks<-sort(unique(r$time));stopifnot(identical(ticks,c(0L,12L,24L,36L)) || all(ticks==c(0,12,24,36)))
 for(t in sort(unique(c(ticks,head(ticks,-1)+diff(ticks)/2))))ln(rep(xp(t),2),c(B,T),'#E6E6E6',.35)
 for(v in seq(0,1,.125))ln(c(L,R),rep(yp(v),2),'#E6E6E6',.35)
 pushViewport(viewport(x=unit(L,'mm'),y=unit(B,'mm'),width=unit(R-L,'mm'),height=unit(T-B,'mm'),
   just=c('left','bottom'),xscale=xr,yscale=yr,clip='on'))
 axes<-measure(paste0(id,'_survival'))
 step<-function(t,v)list(x=c(0,rep(t,each=2)),y=c(1,as.vector(rbind(c(1,head(v,-1)),v))))
 for(g in c('High','Low')){
  a<-z[z$group==g,];a<-a[order(a$time),]
  for(bound in c('lower','upper')){
   v<-a[[bound]];v[is.na(v)&a$survival==0]<-0;stopifnot(!anyNA(v))
   ss<-step(a$time,v)
   grid.lines(unit(ss$x,'native'),unit(ss$y,'native'),gp=gpar(col=colors[[g]],lwd=.6,lty='44'))
  }
 }
 for(g in c('High','Low')){
  a<-z[z$group==g,];a<-a[order(a$time),];ss<-step(a$time,a$survival)
  grid.lines(unit(ss$x,'native'),unit(ss$y,'native'),gp=gpar(col=colors[[g]],lwd=1))
  cens<-a[a$n_censor>0,]
  grid.points(unit(cens$time,'native'),unit(cens$survival,'native'),pch=3,size=unit(1,'mm'),gp=gpar(col=colors[[g]],lwd=.6))
 }
 popViewport();box(L,R,B,T)
 for(v in seq(0,1,.25)){ln(c(L-.5,L),rep(yp(v),2),lwd=.4);txt(sprintf('%.2f',v),L-1,yp(v),6,just='right')}
 for(t in ticks){ln(rep(xp(t),2),c(B-.6,B),lwd=.4);txt(as.character(t),xp(t),26.1,6,just='centre')}
 txt('Event-free survival',2.2,(B+T)/2,6.5,just='centre',rot=90)
 txt(paste0('Log-rank P = ',fmtp(s$logrank_p)),L,22.5,6.3)
 txt(sprintf('HR = %.3f; Cox P = %s',s$HR,fmtp(s$cox_p)),L,19,6.3)
 txt('Number at risk',L,15.9,6.3)
 for(t in head(ticks,-1)+diff(ticks)/2)ln(rep(xp(t),2),c(5.5,14),'#E6E6E6',.35)
 ln(c(L,R),c(9.75,9.75),'#E6E6E6',.35);box(L,R,5.5,14)
 for(g in c('High','Low')){
  yy<-if(g=='High')11.8 else 7.8;txt(g,L-1,yy,6.3,just='right',col=colors[[g]])
  a<-r[r$group==g,]
  for(i in seq_len(nrow(a)))txt(a$n_risk[i],xp(a$time[i]),yy,6.3,just='centre',col=colors[[g]])
 }
 for(t in ticks){ln(rep(xp(t),2),c(4.9,5.5),lwd=.4);txt(as.character(t),xp(t),4,6,just='centre')}
 txt('Time (months)',(L+R)/2,1.6,6.3,just='centre')
 pushViewport(viewport(x=unit(L,'mm'),y=unit(5.5,'mm'),width=unit(R-L,'mm'),height=unit(8.5,'mm'),just=c('left','bottom')))
 rt<-measure(paste0(id,'_risk'));popViewport()
 list(axes,rt)
}
for(marker in names(labels)){
 name<-paste0('Li_CART_EFS_',marker)
 export(name,60,72,function(audit){
  panels<-draw_km(marker)
  if(audit)layout_record(name,60,72,panels,list(list(id='shared_time_axis',panels=c(paste0(marker,'_survival'),paste0(marker,'_risk')))))
 })
}
# True vector comparison sheet; no raster re-import or rescaling of text.
export('Li_CART_four_KMs',124,148,function(audit){
 panels<-list();cols<-list();rows<-list()
 for(i in seq_along(labels)){
  marker<-names(labels)[i];xx<-if(i%%2==1)0 else 64;yy<-if(i<=2)76 else 0
  pushViewport(viewport(x=unit(xx,'mm'),y=unit(yy,'mm'),width=unit(60,'mm'),height=unit(72,'mm'),just=c('left','bottom')))
  panels<-c(panels,draw_km(marker));popViewport()
  cols[[length(cols)+1]]<-list(id=marker,panels=c(paste0(marker,'_survival'),paste0(marker,'_risk')))
 }
 for(r in 1:2)for(kind in c('survival','risk'))rows[[length(rows)+1]]<-list(id=paste0('row',r,kind),panels=paste0(names(labels)[(2*r-1):(2*r)],'_',kind))
 if(audit)layout_record('Li_CART_four_KMs',124,148,panels,cols,rows)
})

# Natural-log effect axes. Drawing on the log scale transforms both limits,
# points and the pooled diamond together; square AREA is proportional to weight.
forest_geometry <- list()
draw_forest <- function(pool_name, audit = FALSE) {
  is_response <- pool_name == 'durability_primary'
  roster <- if (is_response) c('moffitt','zuma7','stahl','li_cart') else c('monti','rchop','chop')
  z <- stats[match(roster,stats$cohort),]
  p <- pools[pools$pool == pool_name,]
  ww <- weights[weights$pool == pool_name,]
  z$weight_percent <- ww$weight_percent[match(z$cohort,ww$cohort)]
  stopifnot(nrow(p)==1L,all(is.finite(z$weight_percent)),abs(sum(z$weight_percent)-100)<1e-8,
            all(abs(log(z$ratio)-z$beta)<1e-10),all(z$CI_low<z$CI_high))
  labels <- c(moffitt='Axi-cel',zuma7='ZUMA-7',stahl='Cologne',li_cart='LymphoMAP',
              monti='GSE34171',rchop='GSE10846-RCHOP',chop='GSE10846-CHOP')
  # The approved forest footprint is 230 x 104 pt, with six-point labels.
  ptmm <- 25.4/72
  text_pt <- function(s,x,y,just='left',bold=FALSE)
    txt(s,x*ptmm,(104-y+2)*ptmm,6,bold,just)
  line_pt <- function(x,y,lwd=.5,lty=1)
    ln(x*ptmm,(104-y)*ptmm,lwd=lwd,lty=lty)
  L <- 89; R <- 163
  limits <- if (is_response) c(-1.5,3) else c(-.2,1.2)
  xp <- function(v) L+(v-limits[1])/diff(limits)*(R-L)
  fmt <- function(v) sprintf('%.2f',v)
  text_pt('Study',3,7,bold=TRUE)
  text_pt(if(is_response)'Log odds ratio [95% CI]' else 'Log hazard ratio [95% CI]',
          227,7,'right',TRUE)
  line_pt(c(1,227),c(13,13))
  line_pt(rep(xp(0),2),c(13,61),lwd=.5,lty='13')
  geometry <- list()
  ys <- if (is_response) 13+8*seq_len(nrow(z)) else c(23,34,45)
  for(i in seq_len(nrow(z))) {
    a <- z[i,]; y <- ys[i]
    lo <- log(a$CI_low); hi <- log(a$CI_high)
    stopifnot(lo>=limits[1],hi<=limits[2])
    side <- sqrt(a$weight_percent)*.38
    line_pt(xp(c(lo,hi)),c(y,y))
    for(v in c(lo,hi))line_pt(rep(xp(v),2),y+c(-.7,.7))
    grid.rect(unit(xp(a$beta)*ptmm,'mm'),unit((104-y)*ptmm,'mm'),
      unit(side*ptmm,'mm'),unit(side*ptmm,'mm'),gp=gpar(fill='black',col=NA))
    text_pt(labels[[a$cohort]],3,y+2)
    text_pt(paste0(fmt(a$beta),' [',fmt(lo),', ',fmt(hi),']'),227,y+2,'right')
    geometry[[i]] <- data.frame(pool=pool_name,cohort=a$cohort,beta=a$beta,
      CI_log_low=lo,CI_log_high=hi,weight_percent=a$weight_percent,
      point_x_pt=xp(a$beta),CI_left_pt=xp(lo),CI_right_pt=xp(hi),
      y_pt=y,square_side_pt=side)
  }
  line_pt(c(1,227),c(51,51))
  lo <- log(p$CI_low); hi <- log(p$CI_high)
  stopifnot(lo>=limits[1],hi<=limits[2])
  xx <- xp(c(lo,p$beta,hi,p$beta)); yy <- 57+c(0,-1.2,0,1.2)
  grid.polygon(unit(xx*ptmm,'mm'),unit((104-yy)*ptmm,'mm'),gp=gpar(fill='black',col=NA))
  text_pt('RE: pooled',3,59)
  text_pt(paste0(fmt(p$beta),' [',fmt(lo),', ',fmt(hi),']'),227,59,'right')
  pv <- if (is_response) sprintf('%.4f',p$p) else formatC(p$p,format='e',digits=2)
  text_pt(paste0('P = ',pv,'; I\u00b2 = ',round(p$I2),'%'),3,66)
  line_pt(c(L,R),c(71,71))
  ticks <- if (is_response) -1:3 else seq(-.2,1.2,.2)
  for(v in ticks) {
    line_pt(rep(xp(v),2),c(71,73))
    text_pt(if(abs(v)<1e-10)'0' else format(round(v,1),trim=TRUE),xp(v),79,'centre')
  }
  text_pt('Better',L-10,79,'right');text_pt('Worse',R+10,79)
  text_pt(if(is_response)'Log odds ratio per IQR of CCL8+CCL13+' else
    'Log hazard ratio per IQR of CCL8+CCL13+',115,87,'centre')
  text_pt(if(is_response)'myeloid cell fraction for non-durable outcome' else
    'myeloid cell fraction for overall survival',115,94.5,'centre')
  text_pt(if(is_response)'6-month or study-defined durable response (DR)' else
    'Fraction of the total deconvolved mixture',115,102,'centre')
  if(audit) {
    forest_geometry[[pool_name]] <<- do.call(rbind,geometry)
    write_json(list(pool=pool_name,log_axis_limits=limits,null=0,better='left',worse='right',
      diamond_log_effects=c(lo,p$beta,hi,p$beta),diamond_x_pt=xx,diamond_y_pt=yy,
      pooled_p=p$p,n=p$n,I2=p$I2,weight_encoding='square area proportional to REML weight'),
      file.path(qa,paste0(pool_name,'_forest_geometry.json')),auto_unbox=TRUE,pretty=TRUE,digits=16)
  }
}
export('Figure5e_chemotherapy_OS',230/72*25.4,104/72*25.4,function(audit)
  draw_forest('chemo_OS',audit))
export('Figure5f_CAR_T_durability',230/72*25.4,104/72*25.4,function(audit)
  draw_forest('durability_primary',audit))
write.csv(do.call(rbind,forest_geometry),file.path(qa,'forest_geometry.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(qa,'R_sessionInfo.txt'))
write_json(list(status='RENDERED',numerical_source='Newly generated analysis output tables',
  KM_dimensions_mm=c(60,72),KM_sheet_dimensions_mm=c(124,148),forest_dimensions_pt=c(230,104),
  minimum_font_pt=6,PNG_dpi=600,models_refitted_in_renderer=FALSE,
  Cologne_reference='oversampled two-cell AND reference; exploratory',
  PDF_text='Editable Helvetica; base PDF avoids a platform-specific Cairo dependency',
  visual_QA='Review PDF/SVG outputs and audit records at final physical size'),
  file.path(qa,'render_record.json'),auto_unbox=TRUE,pretty=TRUE)
