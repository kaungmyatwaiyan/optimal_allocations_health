rm(list=ls())

#############################################################################################################################
### USER INPUTS ###

start_time <- Sys.time() #mark start time

nreps<-1000  #this is the number of replications performed for strategies 4 and 5 (not relevant for strategies 0-3)

g<-0.5 #this is the weight applied to the value component of the objective function for strategies 4 and 5 
#^(see code section 1.2.5 and definitions for strategies 4 and 5 in the paper, e.g. eqn. 11)

#############################################################################################################################
# 1. SETUP

# 1.1 PACKAGES

require(readxl)
require(httr)
require(tidyverse)
require(parallel)
require(packcircles)
require(gridExtra)

# 1.2 FUNCTIONS

# 1.2.1 function to capitalise first letter of each word
simpleCap<-function(x) {
  s<-strsplit(x," ")[[1]]
  paste(toupper(substring(s,1,1)),substring(s,2),
        sep="",collapse=" ")
}

# 1.2.2 function to sample one value at random from vector (including vector of length one)
sample1<-function(x) {
  if (length(x)<=1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}

# 1.2.3 function to perform allocations using either greedy method or random
alloc_alg<-function(mlog,seed,strat) {
  set.seed(seed) #this is because there are random numbers used - it enables the same result to be derived each time (which can help with continuity / debugging)
  acost<-0 #this is the 'accumulated cost', which starts at zero and is added to by the cost of each assigned activity, at each iteration in the below 'while' loop
  while (any(mlog$status==0)) {     #if all of the mlog statuses =1 then all possible activities have been assigned, i.e. the assignment loop needs to terminate
    if (strat %in% 0:3) {    # baseline or strategies 1-3
      mlog.row<-mlog %>%
        filter(elig==1 & status==0 & cost<=tcost-acost)    #this filters rows of the mlog file representing activities that are eligible (see code section 2.3) and
      #have not yet been assigned (i.e. status=0) and whose assignment would not result in the cost constraint being breached
      if (nrow(mlog.row)==0) break   # if there is no activity that meets the condition then terminate the assignment loop
      mlog.row <- mlog.row %>%    # otherwise sample...
        slice(sample1(which(value==max(value))))  #only 1 activity can be selected - this ensures that with the highest 'value' is selected
    } else if (strat %in% 4:5) {    #strategies 4 or 5
      mlog.row<-mlog %>%
        filter(elig==1 & status==0 & cost<=tcost-acost)
      if (nrow(mlog.row)==0) break
      mlog.row <- mlog.row %>%
        sample_n(1)
    }
    if (nrow(mlog.row)==0) break    # if, after applying the above filters, there is no activities that can be assigned then the assignment loop needs to terminate
    mlog$status[which(mlog$pat.ref==mlog.row$pat.ref & mlog$act.ref==mlog.row$act.ref)]<-1   # otherwise, if there is an activity that can be assigned, then flip the status from 0 to 1 
    # to indicate it is assigned (and therefore cannot be considered for future assignments, as above)
    acost<-acost+mlog.row$cost   # update the accumulated cost by the cost of the newly assigned activity
  }
  # when the assignment loop has terminated then add on the seed to the mlog table
  mlog<-mlog %>%
    mutate(seed=seed)
  return(mlog)
}

# 1.2.4 function for selecting best allocations for strategies 4 and 5
alloc_sel_s45<-function(res,strat) {
  tres_a<-res %>%
    group_by(seed) %>%
    summarise(objfn_t1a=sum(status*value),cost=sum(status*cost)) %>%
    mutate(objfn_t1b=(max(objfn_t1a)-objfn_t1a)/(max(objfn_t1a)-min(objfn_t1a))) %>%
    mutate(objfn_t1b=ifelse(is.na(objfn_t1b),0.5,objfn_t1b))
  tres_b<-res %>%
    group_by(seed,pat.ref) %>%
    summarise(value=sum(status*var.opt)) %>%
    summarise(objfn_t2a=var(value)) %>%
    mutate(objfn_t2b=(objfn_t2a-min(objfn_t2a))/(max(objfn_t2a)-min(objfn_t2a))) %>%
    mutate(objfn_t2b=ifelse(is.na(objfn_t2b),0.5,objfn_t2b))
  tres<-left_join(tres_a,tres_b,by="seed") %>%
    mutate(objfn_a=g*objfn_t1a-(1-g)*objfn_t2a) %>%
    mutate(objfn_b=g*objfn_t1b+(1-g)*objfn_t2b)
  tres_out<-res %>%
    #filter(seed==tres$seed[which.max(tres$objfn_a)])
    filter(seed==tres$seed[which.min(tres$objfn_b)])
  return(tres_out)
}

# 1.2.5 function for calculating objective function/objfn for simulated annealing
objectiveFunction <- function(res){
  res_a <- res %>%
    mutate(objfn_t1 = sum(status*value)) #calculate first component of objfn
  res_b <- res %>%
    group_by(pat.ref) %>%
    summarise(ind_value = sum(status*var.opt)) %>%
    mutate(objfn_t2 = var(ind_value)) #calculate second component of objfn
  res<-left_join(res_a,res_b,by='pat.ref') %>%
    mutate(objfn = g*objfn_t1 - (1-g)*objfn_t2) #calculate obj function value
  return(list(res$objfn[1], res_a$objfn_t1[1], res_b$objfn_t2[1])) #returns objfn value, short and long term value, variance)
}

# 1.2.6 Simulated Annealing function for Strategies 4 & 5
simulatedAnnealing<-function(res, seed, niter, step){ 
  '
  Parameters
  
  res - dataframe, which is the best allocation from 1000 replications
  seed - random number generator (RNG) state for random number generation
  niter - number of iterations 
  step - size of the step in reduction of temperature
  '
  set.seed(seed)
  current_res <- res #assign 'res' as current allocation
  best_res <- current_res #current allocation is the best at the moment
  current_fv <- objectiveFunction(current_res)[[1]] #calculate function value of current allocation
  best_fv <- current_fv #current allocation is the best at the moment
  for (k in 1:as.numeric(niter)) {
    Temp <- (1-as.numeric(step))^k #calculate the temperature from k and step
    
    nb_res <- current_res #neighbour allocation will be deduced from current one
    nb_tcost <- sum(nb_res$status*nb_res$cost) #total cost of neighbour/current allocation
    
    nb_alloc <- nb_res %>%
      filter(status==1) #filtering activities included in current allocation
    nb_unalloc <- nb_res %>%
      filter(status==0) #filtering activities NOT included in current allocation
    
    nb_res.A <- nb_alloc %>% 
      sample_n(1) #randomly select one from allocated activities
    
    j = 0
    while (j < 10000){
      nb_res.U <- nb_unalloc %>% 
        sample_n(1) #randomly select one from allocated activities
      if ((nb_res.U$cost + nb_tcost - nb_res.A$cost) <= tcost) {
        #if swapping the allocated and unallocated ones does not violate the cost constraint, the activity will be selected.
        #unless, it will loop until such activity is selected
        break
      }
      j = j + 1
    }
    
    #updating neighbour allocation
    nb_res$status[which(nb_res$pat.ref==nb_res.A$pat.ref & nb_res$act.ref==nb_res.A$act.ref)]<-0
    nb_res$status[which(nb_res$pat.ref==nb_res.U$pat.ref & nb_res$act.ref==nb_res.U$act.ref)]<-1
    
    nb_fv <- objectiveFunction(nb_res)[[1]] #calculate function value of neighbour allocation
    
    if (current_fv <= nb_fv){ #allocation with better fv is selected
      current_res <- nb_res
      current_fv <- nb_fv
      if (best_fv <= current_fv){ #comparison with best allocation
        best_res <- current_res
        best_fv <- current_fv
      }
    } else if (exp((current_fv-nb_fv)/Temp) > runif(1,0,1)){ 
      #if neighbour is not better than current,accept the move with a probability based on the current temperature
      current_res <- nb_res
      current_fv <- nb_fv
    }
  }
  return(list(best_res, best_fv))
}

#############################################################################################################################
# 2. READ DATA AND PRE-PROCESSING

# 2.1 read data
input_url<- "https://raw.githubusercontent.com/kaungmyatwaiyan/optimal_allocations_health/main/inputs_core_seg.xlsx"

temp_file <- tempfile()
GET(input_url, write_disk(temp_file))

x1<-read_excel(temp_file,sheet="attributes")
x2<-read_excel(temp_file,sheet="activity")
x3<-read_excel(temp_file,sheet="cost")

# 2.2 cost
tcost<-x3$total.cost     # the total cost that cannot be exceeded under and circumstances (this would be an 'infeasible solution')

# 2.3 set up raw attributes and activity log - this is never changed after being defined here
aalog_raw<-left_join(x2,x1,by="cohort") %>%
  uncount(volume) %>%     # this duplicates the activity rows by the provided 'volumes' (e.g. if 2x GP consultation demanded, then they'll be two rows)
  mutate(act.ref=1:nrow(.)) %>%
  mutate(val_s0=0,cost_s0=stcp,     # see the paper for these, they should be the same equations (all should be described in the paper)
         val_s1=sto/stcp,cost_s1=stcp,
         val_s2=(sto+lto)/(stcp-ltcp),cost_s2=stcp,
         val_s3=(sto+lto)/(stcp+stcc-ltcp),cost_s3=stcp+stcc,
         val_s4=(sto+lto)/(stcp+stcc-ltcp),cost_s4=stcp+stcc,
         val_s5=(sto+lto)/(stcp+stcc-ltcp),cost_s5=stcp+stcc) %>%
  mutate(elig_s0=ifelse(stcc==0,1,0),
         elig_s1=ifelse(stcc==0,1,0),
         elig_s2=ifelse(stcc==0,1,0),
         elig_s3=1,
         elig_s4=1,
         elig_s5=1)

# 2.4 set up raw main log
mlog_raw_assist<-aalog_raw %>%   # this object simply includes the cumulative count of patients in each of the cohort...
  group_by(cohort) %>%         # ... it is used just to help in the construction of mlog_raw (below)
  summarise(number=unique(number)) %>%
  mutate(cumnumber=cumsum(number)) %>%
  select(cohort,cumnumber)
mlog_raw<-left_join(aalog_raw,mlog_raw_assist,by="cohort") %>%
  select(cohort,act.ref,number,cumnumber) %>%
  uncount(number,.remove=FALSE) %>%
  group_by(cohort,act.ref) %>%
  mutate(pat.ref={cumnumber-number+1}:{cumnumber}) %>%
  ungroup() %>%
  select(pat.ref,act.ref) %>%
  mutate(status=0)
# mlog_raw is a table with as many rows as there are distinct patient activities demanded/needed
# it has 3 columns including the pat.ref (1 to the total number of patients, i.e. the sum of column B on the 'attributes' tab of the excel input file)
# ... the act.ref is a reference to the distinct activities included in aalog_raw (column act.ref)
# ... status is whether the activity is currently assigned or not - the default is zero, and it flips to 1 when assigned in, e.g., alloc_alg in code section 1.2.4

#############################################################################################################################
# 3. PERFORM ALLOCATIONS

strategy_names<-c("baseline","strategy 1","strategy 2","strategy 3","strategy 4","strategy 5")

# 3.1 baseline
aalog_s0<-aalog_raw %>%    # this is taking the raw activity log and selecting only the necessary variables for making assignments
  select(act.ref,val_s0,cost_s0,elig_s0) %>%    # the necessary variables are the value and cost of the activity, and whether it is 'eligible' 
  #(i.e. whether or not it can be considered for assigment)
  rename(value=val_s0,cost=cost_s0,elig=elig_s0)    # standardises column names so it can be used by the generic assignment function
mlog_s0<-left_join(mlog_raw,aalog_s0,by="act.ref")   # joins the above to the mlog_raw, thus creating a table with all possible demanded/needed activities,...
# .... and their various values, costs, and eligibility
res_s0<-alloc_alg(mlog=mlog_s0,seed=1,strat=0)      # this is then passed to the allocation algorithm

# 3.2 strategy 1: greedy algorithm to assign highest SHORT-TERM value activity until cost reached 
aalog_s1<-aalog_raw %>%        #same things happening as above for 3.1
  select(act.ref,val_s1,cost_s1,elig_s1) %>%
  rename(value=val_s1,cost=cost_s1,elig=elig_s1)
mlog_s1<-left_join(mlog_raw,aalog_s1,by="act.ref")
res_s1<-alloc_alg(mlog=mlog_s1,seed=1,strat=1)

# 3.3 strategy 2: greedy algorithm to assign highest SHORT AND LONG-TERM value activity until cost reached
aalog_s2<-aalog_raw %>%
  select(act.ref,val_s2,cost_s2,elig_s2) %>%
  rename(value=val_s2,cost=cost_s2,elig=elig_s2)
mlog_s2<-left_join(mlog_raw,aalog_s2,by="act.ref")
res_s2<-alloc_alg(mlog=mlog_s2,seed=1,strat=2)

# 3.4 strategy 3: greedy algorithm to assign highest SHORT AND LONG TERM value activity APPRECIATING COST-COERSION until cost reached
aalog_s3<-aalog_raw %>%
  select(act.ref,val_s3,cost_s3,elig_s3) %>%
  rename(value=val_s3,cost=cost_s3,elig=elig_s3)
mlog_s3<-left_join(mlog_raw,aalog_s3,by="act.ref")
res_s3<-alloc_alg(mlog=mlog_s3,seed=1,strat=3)

# 3.5 strategy 4: multiple runs of random greedy algorithms
aalog_s4<-aalog_raw %>%
  mutate(var.opt=sto/(stcp+stcc)) %>%
  select(act.ref,val_s4,cost_s4,elig_s4,var.opt) %>%
  rename(value=val_s4,cost=cost_s4,elig=elig_s4)
mlog_s4<-left_join(mlog_raw,aalog_s4,by="act.ref")
cl<-makeCluster(detectCores()-1)
clusterExport(cl=cl,varlist=c("tcost","sample1"),envir=environment())
clusterEvalQ(cl=cl,c(library(dplyr)))
tres_s4<-do.call("bind_rows",parLapply(cl,1:nreps,alloc_alg,mlog=mlog_s4,strat=4))
stopCluster(cl)
alloc_s4<-alloc_sel_s45(res=tres_s4,strat=4)
res_s4_R<-alloc_s4

# 3.5.1 strategy 4: simulated annealing
best_alloc_s4 <- simulatedAnnealing(res_s4_R, 16, 1000, 0.1)
res_s4 <- best_alloc_s4[[1]] #dataframe

# 3.6 strategy 5: multiple runs of random greedy algorithms
aalog_s5<-aalog_raw %>%
  mutate(var.opt=sto) %>%
  select(act.ref,val_s5,cost_s5,elig_s5,var.opt) %>%
  rename(value=val_s5,cost=cost_s5,elig=elig_s5)
mlog_s5<-left_join(mlog_raw,aalog_s5,by="act.ref")
cl<-makeCluster(detectCores()-1)
clusterExport(cl=cl,varlist=c("tcost","sample1"),envir=environment())
clusterEvalQ(cl=cl,c(library(dplyr)))
tres_s5<-do.call("bind_rows",parLapply(cl,1:nreps,alloc_alg,mlog=mlog_s5,strat=5))
stopCluster(cl)
alloc_s5<-alloc_sel_s45(res=tres_s5,strat=5)
res_s5_R<-alloc_s5

# 3.6.1 strategy 5: simulated annealing
best_alloc_s5 <- simulatedAnnealing(res_s5_R, 16, 1000, 0.1)
res_s5 <- best_alloc_s5[[1]] #dataframe

#############################################################################################################################
# 4. CALCULATION DATASETS
# essentially, this code is pulling together the results as generated in code section 3, into a common 'calculation dataset'

res_raw<-do.call("bind_rows",lapply(0:5,function(i) {
  y1<-get(paste0("res_s",i)) %>%
    select(pat.ref,act.ref,status,seed,elig)
  y2<-left_join(y1,aalog_raw %>% 
      select(act.ref,pod1,pod2,name,sto,stcp,stcc,lto,ltcp),by="act.ref") %>%
    select(-act.ref) %>%
    mutate(stv=sto/(stcp+stcc),sltv=(sto+lto)/(stcp+stcc-ltcp),
           expr.demand=ifelse(stcc==0,1,0)) %>%
    mutate(strategy=strategy_names[i+1])
  return(y2)
}))

resA_names<-data.frame(measure_id=paste0("v",1:9),
                       measure=c("Expressed_demand_allocated",
                                 "Unexpressed_demand_allocated",
                                 "All_demand_allocated",
                                 "Short_term_cost_incurred",
                                 "Long_term_cost_avoided",
                                 "Short_term_outcome_improvement",
                                 "Long_term_outcome_deterioration_avoided",
                                 "Short_term_value",
                                 "Short_and_long_term_value")) %>%
  mutate(measure=gsub("_"," ",measure)) %>%
  mutate(measure=factor(measure,levels=measure))

resA_values<-res_raw %>%
  group_by(strategy,name,pod1,pod2,pat.ref) %>%
  summarise(v1=sum(status==1 & expr.demand==1),
            v2=sum(status==1 & expr.demand==0),
            v3=sum(status==1),
            v4=sum(stcp*status+stcc*status),
            v5=sum(ltcp*status),
            v6=sum(sto*status),v7=sum(lto*status),
            v8=sum(stv*status),v9=sum(sltv*status),
            nreps=length(unique(seed))) %>%
  mutate_at(vars(starts_with("v")),~./nreps) %>%
  select(-nreps) %>%
  pivot_longer(cols=-c(strategy,name,pod1,pod2,pat.ref),names_to="measure_id",
               values_to="values")

resA_poss1<-res_raw %>%
  group_by(strategy,name,pod1,pod2,pat.ref) %>%
  summarise(v1=sum(expr.demand==1),v2=sum(expr.demand==0),v3=sum(expr.demand>=0),
            v4=sum(stcp+stcc),v5=sum(ltcp),
            v6=sum(sto),v7=sum(lto),
            v8=sum(stv),v9=sum(sltv),nreps=length(unique(seed))) %>%
  mutate_at(vars(starts_with("v")),~./nreps) %>%
  select(-nreps) %>%
  pivot_longer(cols=-c(strategy,name,pod1,pod2,pat.ref),
               names_to="measure_id",values_to="poss1")

resA_poss2<-res_raw %>%
  group_by(strategy,name,pod1,pod2,pat.ref) %>%
  filter(elig==1) %>%
  summarise(v1=sum(expr.demand==1),v2=sum(expr.demand==0),v3=sum(expr.demand>=0),
            v4=sum(stcp+stcc),v5=sum(ltcp),
            v6=sum(sto),v7=sum(lto),
            v8=sum(stv),v9=sum(sltv),nreps=length(unique(seed))) %>%
  mutate_at(vars(starts_with("v")),~./nreps) %>%
  select(-nreps) %>%
  pivot_longer(cols=-c(strategy,name,pod1,pod2,pat.ref),
               names_to="measure_id",values_to="poss2")

resA<-left_join(resA_values,resA_poss1,
                by=c("strategy","name","pod1","pod2","pat.ref","measure_id"))%>%
  left_join(resA_poss2,
            by=c("strategy","name","pod1","pod2","pat.ref","measure_id")) %>%
  left_join(resA_names,by="measure_id") %>%
  select(-measure_id) %>%
  mutate(poss1=replace_na(poss1,0)) %>%
  mutate(poss2=replace_na(poss2,0))

# calculating short and long term value and variance of allocations by the former method vs simulated annealing
val_s4_SA <- objectiveFunction(res_s4)[[2]]
var_s4_SA <- objectiveFunction(res_s4)[[3]]
val_s4_R <- objectiveFunction(res_s4_R)[[2]]
var_s4_R <- objectiveFunction(res_s4_R)[[3]]

val_s5_SA <- objectiveFunction(res_s5)[[2]]
var_s5_SA <- objectiveFunction(res_s5)[[3]]
val_s5_R <- objectiveFunction(res_s5_R)[[2]]
var_s5_R <- objectiveFunction(res_s5_R)[[3]]

#############################################################################################################################
# 5. SUMMARISE RESULTS
# producing the plots etc from the calculation datasets generated in code section 4

######################################################
# plot 2: cohort and activity level

pdf("plot_2.0_allocations_by_strategy.pdf",height=6,width=8)
td1<- subset(resA, strategy=='baseline' | 
               strategy == 'strategy 4' | 
               strategy == 'strategy 5') %>%
  group_by(strategy,measure,pod1) %>%
  summarise(values=sum(values),poss1=sum(poss1),poss2=sum(poss2))
td2<-subset(resA, strategy=='baseline' | 
              strategy == 'strategy 4' | 
              strategy == 'strategy 5') %>%
  group_by(strategy,measure) %>%
  summarise(values=sum(values),poss1=sum(poss1),poss2=sum(poss2))
print(td1 %>%
        ggplot(aes(x=strategy,y=values)) +
        geom_bar(stat="identity",aes(fill=pod1)) +
        facet_wrap(~measure,scales="free_y",nrow=4) +
        labs(title="Allocations") +
        xlab("Strategy") +
        theme(legend.position="bottom",
              legend.title=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_text(angle=25,hjust=1),
              strip.text.x=element_text(size=8))
)
dev.off()

######################################################
# plot 3: variances (person-level)

pdf("plot_3.0_variance_by_strategy.pdf",height=6,width=8)
print(subset(resA, strategy=='baseline' | 
               strategy == 'strategy 4' | 
               strategy == 'strategy 5') %>%
        group_by(strategy,pat.ref,measure) %>%
        summarise(values=sum(values)) %>%
        group_by(strategy,measure) %>%
        summarise(values=var(values)) %>%
        ggplot(aes(x=strategy,y=values)) +
        geom_bar(stat="identity") +
        facet_wrap(~measure,scales="free_y",nrow=4) +
        labs(title="Variance between individuals") +
        xlab("Strategy") +
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_text(angle=25,hjust=1),
              strip.text.x=element_text(size=8))
)
dev.off()

end_time <- Sys.time() #mark end time
end_time - start_time #time spent


