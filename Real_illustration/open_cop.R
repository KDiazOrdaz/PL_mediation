cop_df <- read.csv('Datasets/copers.csv')

head(cop_df)
names(cop_df)

cop_treat <- subset(cop_df,treat==1)[,-c(1,18)]
dim(cop_treat)

table(cop_treat$sessions_attended)
require(ggplot2)

ggplot() + geom_bar(data=cop_treat,aes(x=sessions_attended))

cop_treat$X <- as.numeric(cop_treat$sessions_attended >= 12)
cop_treat$M <- cop_treat$PSEQ_overall_12wk
cop_treat$Y <- cop_treat$CPG_overall_12m

table(cop_treat$X)

Z_formula = "~Gender+age+hads_dep_base+hads_anx_base+PSEQ_overall_base+
HEIQ_base+CPAQ_base+Housing+English+Employment+
Education+employed_ftedu+centre"

Z_formula = "~Gender+poly(age,degree=2)+hads_anx_base+employed_ftedu*centre+hads_dep_base+
    PSEQ_overall_base+
    HEIQ_base+CPAQ_base"

mod1 = glm(paste0('X',Z_formula),family=binomial,data=cop_treat)
summary(mod1)

cop_treat$p_x = predict(mod1,cop_treat,type = "response")

ggplot()+geom_density(aes(x=p_x,fill=as.character(X)),data=cop_treat,alpha=0.4)+xlim(0,1)



mod2 = glm(paste0('M',Z_formula,'+X'),data=cop_treat)
summary(mod2)
ggplot()+geom_point(aes(x=p_x,y=predict(mod2,cop_treat),
                        color=as.character(X)),data=cop_treat,alpha=0.9)


names(cop_treat)

ind_mis = is.na(cop_treat$p_x)|is.na(cop_treat$M)|is.na(cop_treat$Y)
sum(ind_mis)
sum(is.na(cop_treat$p_x))
sum(is.na(cop_treat$M)|is.na(cop_treat$Y))
sum(ind_mis) - sum(is.na(cop_treat$M)|is.na(cop_treat$Y))
##Missing data is a bit annoying
# of 384 people, 74 are missing either M,Y (M/Y=54) or baseline covariate

missings = cop_treat[ind_mis,]
#View(missings)

cop_treat_complete = cop_treat[!ind_mis,]
dim(cop_treat_complete)




