# Meta analysis
library(meta)
library(metafor)

df <- read.csv("datatable_all.csv")
df$AuthYear <- paste(df$Author,df$Year, sep=", ")
#df$ntot <- df$Pop_Drug + df$Pop_Placebo

# Venn diagram
library(VennDiagram)
A <- subset(df$AuthYear, df$Borst == "yes")
B <- subset(df$AuthYear, df$Xu == "yes")
C <- subset(df$AuthYear, df$Corona == "yes")
D <- subset(df$AuthYear, df$Calof == "yes")
E <- subset(df$AuthYear, df$Haddad == "yes")

venn.plot <- venn.diagram(
    x = list(
        Borst = A,
        Xu = B,
        Corona = C,
        Calof = D,
        Haddad = E
    ),
    filename = "venn5.tiff",
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1","orange"),
    alpha = 0.50,
    label.col = rep("black",31),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4","orange"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 270,
    margin = 0.2
);

### ----------------------------------------------------------------------------
### Borst
### ----------------------------------------------------------------------------

# figures are 12" by 9"

### Borst random-effects conditional Poisson model
res <- rma.glmm(measure = "IRR",
                x1i = Events_Drug,
                t1i = Pop_Drug,
                x2i = Events_Placebo,
                t2i = Pop_Placebo,
                mods = ~Mode - 1,
                data = subset(df, df$Borst == "yes"),
                model = "CM.EL",
                drop00 = F,
                tdist = T,
                slab = AuthYear)

res.all <- rma.glmm(measure = "IRR",
                    x1i = Events_Drug,
                    t1i = Pop_Drug,
                    x2i = Events_Placebo,
                    t2i = Pop_Placebo,
                    data = subset(df,df$Borst == "yes"),
                    model = "CM.EL",
                    drop00 = F,
                    tdist = T,
                    slab = AuthYear)

print(res, digits=3)

forest(res,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim = c(-3,3),
       ylim = c(-10,40),
       cex = .8)
addpoly(c(res$b[c(1,4,2,3)],res.all$b),
        sei = c(res$se[c(1,4,2,3)], res.all$se),
        atransf = exp,
        mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)],
        cex = .8)

### Borst random-effects conditional Poisson model, with 3 papers removed
removed  <- c(5,10,13)
#res <- rma.glmm(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, data=df, model="CM.EL", drop00=F, tdist=T, slab=AuthYear)

res <- rma.glmm(measure="IRR",
                x1i=Events_Drug,
                t1i=Pop_Drug,
                x2i=Events_Placebo,
                t2i=Pop_Placebo,
                mods=~Mode-1,
                data=subset(df[-removed,],df[-removed,]$Borst=="yes"),
                model="CM.EL",
                drop00=F,
                tdist=T,
                slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df[-removed,],df[-removed,]$Borst=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
#res <- rma(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, mods=~Mode-1, data=df, drop00=F, knha=T, slab=AuthYear)
print(res, digits=3)
print(res.all, digits=3)
forest(res,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
addpoly(c(res$b[c(1,4,2,3)],res.all$b),
        sei = c(res$se[c(1,4,2,3)],res.all$se),
        atransf = exp,
        mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)],
        cex=.8)


### ----------------------------------------------------------------------------
### Xu
### ----------------------------------------------------------------------------

# figures are 12" by 9"

### Xu random-effects conditional Poisson model
res <- rma.glmm(measure="IRR",
                x1i=Events_Drug,
                t1i=Pop_Drug,
                x2i=Events_Placebo,
                t2i=Pop_Placebo,
                mods=~Mode-1,
                data=subset(df,df$Xu=="yes"),
                model="CM.EL",
                drop00=F,
                tdist=T,
                slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df,df$Xu=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
print(res, digits=3)
forest(res,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
addpoly(c(res$b[c(1,4,2,3)],res.all$b),
        sei = c(res$se[c(1,4,2,3)],res.all$se),
        atransf = exp,
        mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)],
        cex=.8)

### Xu random-effects conditional Poisson model, with 3 papers removed
removed  <- c(5,10,13)
res <- rma.glmm(measure="IRR",
                x1i=Events_Drug,
                t1i=Pop_Drug,
                x2i=Events_Placebo,
                t2i=Pop_Placebo,
                mods=~Mode-1,
                data=subset(df[-removed,],df[-removed,]$Xu=="yes"),
                model="CM.EL",
                drop00=F,
                tdist=T,
                slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df[-removed,],df[-removed,]$Xu=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
print(res, digits=3)
print(res.all, digits=3)
forest(res,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
addpoly(c(res$b[c(1,4,2,3)],res.all$b),
        sei = c(res$se[c(1,4,2,3)],res.all$se),
        atransf = exp,
        mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)],
        cex=.8)


### ----------------------------------------------------------------------------
### Corona
### ----------------------------------------------------------------------------

# figures are 12" by 9"

### Corona random-effects conditional Poisson model
#res <- rma.glmm(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, mods=~Mode-1, data=subset(df,df$Corona=="yes"), model="CM.EL", drop00=F, tdist=T, slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df,df$Corona=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
#print(res, digits=3)
forest(res.all,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
#addpoly(c(res$b[c(1,4,2,3)],res.all$b), sei = c(res$se[c(1,4,2,3)],res.all$se), atransf = exp, mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)], cex=.8)

### Corona random-effects conditional Poisson model, with 3 papers removed
removed  <- c(5,10,13)
#res <- rma.glmm(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, mods=~Mode-1, data=subset(df[-removed,],df[-removed,]$Corona=="yes"), model="CM.EL", drop00=F, tdist=T, slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df[-removed,],df[-removed,]$Corona=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
print(res, digits=3)
print(res.all, digits=3)
forest(res.all,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
#addpoly(c(res$b[c(1,4,2,3)],res.all$b), sei = c(res$se[c(1,4,2,3)],res.all$se), atransf = exp, mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)], cex=.8)



### ----------------------------------------------------------------------------
### Haddad
### ----------------------------------------------------------------------------

# figures are 12" by 9"

### Haddad random-effects conditional Poisson model
#res <- rma.glmm(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, mods=~Mode-1, data=subset(df,df$Corona=="yes"), model="CM.EL", drop00=F, tdist=T, slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df,df$Haddad=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
#print(res, digits=3)
forest(res.all,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
#addpoly(c(res$b[c(1,4,2,3)],res.all$b), sei = c(res$se[c(1,4,2,3)],res.all$se), atransf = exp, mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)], cex=.8)

### Haddad random-effects conditional Poisson model, with 3 papers removed
removed  <- c(5,10,13)
#res <- rma.glmm(measure="IRR", x1i=Events_Drug, t1i=Pop_Drug, x2i=Events_Placebo, t2i=Pop_Placebo, mods=~Mode-1, data=subset(df[-removed,],df[-removed,]$Corona=="yes"), model="CM.EL", drop00=F, tdist=T, slab=AuthYear)
res.all <- rma.glmm(measure="IRR",
                    x1i=Events_Drug,
                    t1i=Pop_Drug,
                    x2i=Events_Placebo,
                    t2i=Pop_Placebo,
                    data=subset(df[-removed,],df[-removed,]$Haddad=="yes"),
                    model="CM.EL",
                    drop00=F,
                    tdist=T,
                    slab=AuthYear)
print(res, digits=3)
print(res.all, digits=3)
forest(res.all,
       atransf = exp,
       at = log(c(0.05, 0.25, 1, 4, 20)),
       xlim = c(-9, 7),
       alim=c(-3,3),
       ylim=c(-10,40),
       cex=.8)
#addpoly(c(res$b[c(1,4,2,3)],res.all$b), sei = c(res$se[c(1,4,2,3)],res.all$se), atransf = exp, mlab = c("gel","injection","oral","patch","all")[c(1,4,2,3,5)], cex=.8)


