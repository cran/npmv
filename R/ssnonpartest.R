#Function to perform the subset algorithm.  Algorithm is comprised of 3 steps.  Step 1: Global Hypothesis Test.  Step 2: Test all a-1 (or p-1) subsets.  
# Step 3: Test all remaing subsets per closed multiple testing principle


ssnonpartest <- function(formula,data,alpha=.05,test=c(1,0,0,0),factors.and.variables=FALSE){

#Checks to see if formula
  if(!is(formula,"formula")){
    return('Error: Please give a formula')
  }  

#Checks to ensure only one test is requested
  if(sum(test)!=1){
    return('Error:Please specify a single test')
  }
  
#Creates the data frame
  formula=Formula(formula)
  frame=model.frame(formula,data=data)
  
   
#Checks for missing data
   if(sum(is.na(frame))>0)
   {
      return('Error: Missing Data')
   }
  
#Assigns group variable and response variables
  groupvar.location=length(frame[1,])
  groupvar=names(frame)[groupvar.location]
  vars=names(frame)[1:(groupvar.location-1)]

#Changes factor levels to a factor if not already
   if(!is.factor(frame[,groupvar]))
   {
      frame[,groupvar] <- factor(frame[,groupvar])
   }
   
#Gives levels of factor
   levels=levels(frame[,groupvar])

# Orders data by group
   o <- order(frame[,groupvar])
   frame <- frame[o,]

# Compute a-number of factors and --number of variables
   p <- length(vars)
   a <- length(levels(frame[,groupvar]))

#Defines logical variable that tells when to exit
  exit=FALSE
  
# Step 1: Test the Global Hypothesis 
 
   base <- basenonpartest(frame,groupvar,vars,tests=test)
   if (test[1]==1){testpval=base$pvalanova}
   if (test[2]==1){testpval=base$pvalLH}
   if (test[3]==1){testpval=base$pvalBNP}
   if (test[4]==1){testpval=base$pvalWL}
   
   if( testpval < alpha) {cat('The Global Hypothesis is significant, subset algorithm will continue \n')
   } else {return ('The Global Hypothesis is not significant, subset algorithm will not continue')}

#Algorithm for when p>a
if(p>a || factors.and.variables==TRUE){
cat('\n~Performing the Subset Algorithm based on Factor levels~\nThe Hypothesis between factor levels ', levels, 'is significant  \n')
# Step 2: Subset Algorithm for Testing Factor Levels with 1 factor removed, returns matrix where the rows are subsets(size a-1) that are significant
#Creates new matrix of all possible intersections and finds which are signficant, creates new matrix of significant subsets

#Test each of the a-1 subgroups, creates list of a-1 groups which are significant
step2subsets=vector("list",a)
for(i in 1:a)
  {
  subsetframe <-subset(frame, frame[,groupvar] != levels[i])
  subsetframe<- droplevels(subsetframe)
  groupvarsub=names(subsetframe)[groupvar.location]
  base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
   if (test[1]==1){testpval=base$pvalanova}
   if (test[2]==1){testpval=base$pvalLH}
   if (test[3]==1){testpval=base$pvalBNP}
   if (test[4]==1){testpval=base$pvalWL}
   

   if( testpval < alpha) {step2subsets[[i]]=levels(subsetframe[,groupvarsub])}
    else{step2subsets[[i]]=NA}
   if( testpval < alpha) {cat('The Hypothesis between factor levels ', siglevels=levels[-i], 'is significant  \n')}
  }

step2subsets=step2subsets[!is.na(step2subsets)]

if (length(step2subsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
if (length(step2subsets)<= 1 && factors.and.variables==TRUE){
  cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
  exit=TRUE}

#Creates a new list of all the possible intersections of the previous significant subsets
if(exit==FALSE){
step2subsetcount=length(step2subsets)
num.intersections=((step2subsetcount-1)*step2subsetcount)/2
newsubsets=vector("list",num.intersections)
k=1   
   for(i in 1:(step2subsetcount-1)){
     h=i+1
      for(j in h:step2subsetcount){
        newsubsets[[k]]=intersect(step2subsets[[i]],step2subsets[[j]])
        k=k+1
      }
   }  
  

#Creates a list of significant subsets and non-significant subsets
newsubsetcount=length(newsubsets)
sigfactorsubsets=vector("list",newsubsetcount)
nonsigfactorsubsets=vector("list",newsubsetcount)
for(i in 1:newsubsetcount)
   {
     subsetstotest=as.factor(newsubsets[[i]])
     subsetframe <-subset(frame, frame[,groupvar] %in% subsetstotest)
     subsetframe<- droplevels(subsetframe)
     groupvarsub=names(subsetframe)[groupvar.location]
     base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
     if (test[1]==1){testpval=base$pvalanova}
     if (test[2]==1){testpval=base$pvalLH}
     if (test[3]==1){testpval=base$pvalBNP}
     if (test[4]==1){testpval=base$pvalWL}
     
     if( testpval > alpha) {nonsigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{nonsigfactorsubsets[[i]]=NA}
     if( testpval < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else{sigfactorsubsets[[i]]=NA}
     if( testpval < alpha) {cat('The Hypothesis between factor levels ', newsubsets[[i]], 'is significant  \n')}
   }
   nonsigfactorsubsets=nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
   sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]

   if (length(sigfactorsubsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
   if (length(sigfactorsubsets)<= 1 && factors.and.variables==TRUE){
     cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
     exit=TRUE}
}
#Step 3: Reads in the significant and non significant subsets, intersects the significant sets, checks to see if intersection is 
# non-significant matrix, for all intersections that are not in non-significant matrix creates a new matrix of subsets to check. 
#Creates new matrix of significant and non significant sets
   
#Creates a new list of all the possible intersections of the previous significant subsets

number.elements=a-2
for(l in 1:(a-4)) {  
if(exit==FALSE){  
newsubsetcount=length(sigfactorsubsets)
rows=((newsubsetcount-1)*newsubsetcount)/2
newsubsets=vector("list",rows)
k=1   
   for(i in 1:(newsubsetcount-1)){
     h=i+1
     for(j in h:newsubsetcount){
       newsubsets[[k]]=intersect(sigfactorsubsets[[i]],sigfactorsubsets[[j]])
       k=k+1
     }
   } 
newsubsets=unique(newsubsets)

#Checks to see if intsections are in non-significant subsets, assigns NA to subsets that are in non-significant sets
if(length(nonsigfactorsubsets)>0){
   for(i in 1:length(nonsigfactorsubsets))
   {
     for(j in 1:length(newsubsets))
     {
       if(sum(!(newsubsets[[j]] %in% nonsigfactorsubsets[[i]]))==0){newsubsets[[j]]=NA}
     }
   }
}  
#Removes duplicates
newsubsets=unique(newsubsets)

#Removes subsets of size less than current size (1:a-2)
for(i in 1:length(newsubsets))
  {
  if(length(newsubsets[[i]])<(number.elements-1)){newsubsets[[i]]=NA}
  }
# Removes NA
newsubsets=newsubsets[!is.na(newsubsets)]

#Creates a list of significant subsets and non-significant subsets
newsubsetcount=length(newsubsets)
sigfactorsubsets=vector('list',newsubsetcount)
nonsigfactorsubsets=vector('list',newsubsetcount)

   for(i in 1:newsubsetcount)
   {
     subsetstotest=as.factor(newsubsets[[i]])
     subsetframe <-subset(frame, frame[,groupvar] %in% subsetstotest)
     subsetframe<- droplevels(subsetframe)
     groupvarsub=names(subsetframe)[groupvar.location]
     base <- basenonpartest(subsetframe,groupvarsub,vars,tests=test)
     if (test[1]==1){testpval=base$pvalanova}
     if (test[2]==1){testpval=base$pvalLH}
     if (test[3]==1){testpval=base$pvalBNP}
     if (test[4]==1){testpval=base$pvalWL}
     
     if( testpval > alpha) {nonsigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else {nonsigfactorsubsets[[i]]=NA}
     if( testpval < alpha) {sigfactorsubsets[[i]]=levels(subsetframe[,groupvarsub])}else {sigfactorsubsets[[i]]=NA}
     if( testpval < alpha) {cat('The Hypothesis between factor levels ', newsubsets[[i]], 'is significant  \n')}
   }
   
   nonsigfactorsubsets=nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
   sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]
if (length(sigfactorsubsets)<= 1 && factors.and.variables==FALSE){return(cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
if (length(sigfactorsubsets)<= 1 && factors.and.variables==TRUE){
  cat('All appropriate subsets using factor levels have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha,'\n')
  exit==TRUE}
}
number.elements=number.elements-1
}
}
  
if(p<=a ||factors.and.variables==TRUE){
cat('\n~Performing the Subset Algorithm based on Response Variables~ \nThe Hypothesis involving response variables ', vars, 'is significant  \n')
#Algorithm for p<=a
  # Step 2: Subset Algorithm for Testing Different response variabless with 1 varible removed, returns matrix where the rows are subsets(size p-1) that are significant
  #Creates new matrix of all possible intersections and finds which are signficant, creates new matrix of significant subsets
  
  #Test each of the p-1 subgroups, creates list of p-1 groups which are significant
  step2subsets=vector("list",p)
  for(i in 1:p)
  {
    base <- basenonpartest(frame,groupvar,vars[-i],tests=test)
    if (test[1]==1){testpval=base$pvalanova}
    if (test[2]==1){testpval=base$pvalLH}
    if (test[3]==1){testpval=base$pvalBNP}
    if (test[4]==1){testpval=base$pvalWL}
    
    
    if( testpval < alpha) {step2subsets[[i]]=vars[-i]}
    else{step2subsets[[i]]=NA}
    if( testpval < alpha) {cat('The Hypothesis involving response variables ', sigvariables=vars[-i], 'is significant \n')}
  }
  
  step2subsets=step2subsets[!is.na(step2subsets)]
  
  if (length(step2subsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  
  
  #Creates a new list of all the possible intersections of the previous significant subsets
  
  step2subsetcount=length(step2subsets)
  num.intersections=((step2subsetcount-1)*step2subsetcount)/2
  newsubsets=vector("list",num.intersections)
  k=1   
  for(i in 1:(step2subsetcount-1)){
    h=i+1
    for(j in h:step2subsetcount){
      newsubsets[[k]]=intersect(step2subsets[[i]],step2subsets[[j]])
      k=k+1
    }
  }  
  
  
  #Creates a list of significant subsets and non-significant subsets
  newsubsetcount=length(newsubsets)
  sigfactorsubsets=vector("list",newsubsetcount)
  nonsigfactorsubsets=vector("list",newsubsetcount)
  for(i in 1:newsubsetcount)
  {

    base <- basenonpartest(frame,groupvar,vars=newsubsets[[i]],tests=test)
    if (test[1]==1){testpval=base$pvalanova}
    if (test[2]==1){testpval=base$pvalLH}
    if (test[3]==1){testpval=base$pvalBNP}
    if (test[4]==1){testpval=base$pvalWL}
    
    if( testpval > alpha) {nonsigfactorsubsets[[i]]=newsubsets[[i]]}else{nonsigfactorsubsets[[i]]=NA}
    if( testpval < alpha) {sigfactorsubsets[[i]]=newsubsets[[i]]}else{sigfactorsubsets[[i]]=NA}
    if( testpval < alpha) {cat('The Hypothesis involving response variables ', newsubsets[[i]], 'is significant  \n')}
  }
  nonsigfactorsubsets=nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
  sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]
  if (length(sigfactorsubsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  
  #Step 3: Reads in the significant and non significant subsets, intersects the significant sets, checks to see if intersection is 
  # non-significant list, for all intersections that are not in non-significant matrix creates a new matrix of subsets to check. 
  #Creates new list of significant and non significant sets
  
  #Creates a new list of all the possible intersections of the previous significant subsets
  number.elements=p-2
  for(l in 1:(p-4)) {  
    
    newsubsetcount=length(sigfactorsubsets)
    rows=((newsubsetcount-1)*newsubsetcount)/2
    newsubsets=vector("list",rows)
    k=1   
    for(i in 1:(newsubsetcount-1)){
      h=i+1
      for(j in h:newsubsetcount){
        newsubsets[[k]]=intersect(sigfactorsubsets[[i]],sigfactorsubsets[[j]])
        k=k+1
      }
    } 
    newsubsets=unique(newsubsets)
    
    #Checks to see if intsections are in non-significant subsets, assigns NA to subsets that are in non-significant sets
    if(length(nonsigfactorsubsets)>0){
      for(i in 1:length(nonsigfactorsubsets))
      {
        for(j in 1:length(newsubsets))
        {
          if(sum(!(newsubsets[[j]] %in% nonsigfactorsubsets[[i]]))==0){newsubsets[[j]]=NA}
        }
      }
    }  
    #Removes duplicates
    newsubsets=unique(newsubsets)
    
    #Removes subsets of size less than current size (1:p-2)
    for(i in 1:length(newsubsets))
    {
      if(length(newsubsets[[i]])<(number.elements-1)){newsubsets[[i]]=NA}
    }
    # Removes NA
    newsubsets=newsubsets[!is.na(newsubsets)]
    
    #Creates a list of significant subsets and non-significant subsets
    newsubsetcount=length(newsubsets)
    sigfactorsubsets=vector('list',newsubsetcount)
    nonsigfactorsubsets=vector('list',newsubsetcount)
    
    for(i in 1:newsubsetcount)
    {
     
      base <- basenonpartest(frame,groupvar,vars=newsubsets[[i]],tests=test)
      if (test[1]==1){testpval=base$pvalanova}
      if (test[2]==1){testpval=base$pvalLH}
      if (test[3]==1){testpval=base$pvalBNP}
      if (test[4]==1){testpval=base$pvalWL}
      
      if( testpval > alpha) {nonsigfactorsubsets[[i]]=newsubsets[[i]]}else {nonsigfactorsubsets[[i]]=NA}
      if( testpval < alpha) {sigfactorsubsets[[i]]=newsubsets[[i]]}else {sigfactorsubsets[[i]]=NA}
      if( testpval < alpha) {cat('The Hypothesis involving response variables ', newsubsets[[i]], 'is significant  \n')}
    }
    
    nonsigfactorsubsets=nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
    sigfactorsubsets=sigfactorsubsets[!is.na(sigfactorsubsets)]
    number.elements=number.elements-1
    if (length(sigfactorsubsets)<= 1){return(cat('All appropriate subsets using response variables have been checked using a closed multiple testing procedure, which controls the maximum overall type I error rate at alpha=',alpha))}
  } 
  
  
  
}
}   

