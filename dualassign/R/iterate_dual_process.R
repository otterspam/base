#' Assigns a task to a worker based on the dual variables of a linear program
#'
#' @param lpo list object containing a variety of information about dual linear program
#' @return updated list object with the new assignment and cost increment included
#' @export
#' @examples
#' lpo=list(
#'       n=n,
#'       obj=dual.obj,
#'       rhs=dual.rhs,
#'       mat=dual.matrix,
#'       dir=dual.dir,
#'       ass=c(rep(0,n[3])),
#'       bl=c(rep(0,n[1])),
#'       istar=NULL,
#'       jstar=NULL,
#'       cost=0,
#'       exit=0,
#'       count=0)
#' lpo2=iterate_dual_process(lpo)

iterate_dual_process=function(lpo){
  lpo$istar=0
  lpo$jstar=0
  lp.object=lpSolve::lp(dir="max",lpo$obj,lpo$mat,lpo$dir,lpo$rhs)
  # lpo$mat in block form
  #    theta             eta
  # x [ -I  | perform %*% diag(P[1,]) ]
  # x [ -I  | perform %*% diag(P[2,]) ]
  #...
  # x [ -I  | perform %*% diag(P[K,]) ]
  # y [  0  | diag(otperform)         ]
  eta.index=c((lpo$n[3]+1):(lpo$n[3]+lpo$n[1]))
  theta.index=c(1:lpo$n[3])
  x.index=c(1:(lpo$n[3]*lpo$n[2]))
  y.index=c((lpo$n[3]*lpo$n[2]+1):(lpo$n[3]+lpo$n[2]+lpo$n[1]))

  old.assignment=0
  old.demand=lpo$obj[eta.index]

  eta=lp.object$solution[eta.index]
  demand=lpo$obj[eta.index]
  selection=lpo$mat[,eta.index]%*%eta/lpo$rhs
  optass=which.max(selection)
  #deduct work done from work remaining and put it as the demand going forward
  demand=demand-lpo$mat[optass,eta.index]
  lpo$obj[eta.index]=demand
  if (optass > (lpo$n[3]*lpo$n[2])){
    # this indicates a backlog selection, work out which queue is required
    lpo$istar=0
    lpo$jstar=optass-(lpo$n[3]*lpo$n[2])
    # increment the appropriate backlog and deduct the work done from the demand
    lpo$bl[lpo$jstar]=lpo$bl[lpo$jstar]+1
  }  else {
    # this indicates a normal assignment; determine who is assigned where
    lpo$istar=(optass-1)%%(lpo$n[3])+1
    lpo$jstar=(optass-lpo$istar)/lpo$n[3]+1
    old.assignment=lpo$ass[lpo$istar]
    lpo$ass[lpo$istar]=lpo$jstar
    # set performance rates for the assigned person to be zero as they can no longer work on other things
    unassigned=c(rep(lpo$ass==0,lpo$n[2]))
    lpo$mat[x.index,eta.index]=lpo$mat[x.index,eta.index]*matrix(c(rep(unassigned,lpo$n[1])),ncol=lpo$n[1])
  }
  #add cost of this assignment
  lpo$cost=lpo$cost+lpo$rhs[optass]
  #increment count
  lpo$count=lpo$count+1


  #determine whether to continue: exit codes
  #exit code 1: all demand met this is the good one
  if (max(demand) <= 0){
    lpo$exit=1
  }
  #exit code 2: linear program infeasible or unbounded
  if (lp.object$status > 0){
    lpo$exit=2
  }
  #exit code 3: reassigning someone already assigned
  if (old.assignment>0){
    lpo$exit=3
  }
  #exit code 4: assignment did not reduce demand
  if (max(old.demand>demand)<1){
    lpo$exit=4
  }
  lpo
}
