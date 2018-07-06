#' Assignment procedure with overtime/backlog. Utilizes dual LP for approximate sequential assignment.
#'
#' @param ss list containing state space for the assignment problem. Users are given assignments that split their time among tasks according to P.
#' Users have performance metric for each task. otcost is the relative cost of sending a job to backlog/overtime.
#' otperform is the performance metric for this overtime/backlogged task
#' @return state space updated with assignments and backlog
#' @export
#' @examples
#' ss=list(
#'   tasks=c("first","second","third"),
#'   assigns=c("A","B","c"),
#'   users=c(1:20),
#'   perform=matrix(rgamma(60,5,1),ncol=3),
#'   costs=matrix(c(rep(1,60)),ncol=3),
#'   P=matrix(c(1,0.6,0,0,0.4,0.5,0,0,0.5),ncol=3),
#'   demand=c(300,200,150),
#'   otcost=3,
#'   otperform=c(5,5,5),
#'   solution=NULL,
#'   backlog=NULL,
#'   istar=NULL,
#'   jstar=NULL,
#'   status=c("incomplete")
#' )
#'
#' ss2=dual_assignment(ss)
#' ss2$solution
#' ss2$backlog

dual_assignment=function(ss){
  # set up linear programming parts
  n=c(length(ss$tasks),length(ss$assigns),length(ss$users))
  if (is.null(ss$costs)){ss$costs=matrix(c(rep(1,n[3]*n[2])),ncol=n[2])}
  dual.rhs=c(as.vector(ss$costs),rep(ss$otcost,n[1]))
  dual.obj=c(rep(-1,n[3]),ss$demand)
  dual.dir=c(rep("<=",n[3]*n[2]+n[1]))
  # Construct matrix for dual: In block form
  #          theta             eta
  # assign1 [ -I  | perform %*% diag(P[1,]) ]
  # assign2 [ -I  | perform %*% diag(P[2,]) ]
  #...
  # assignK [ -I  | perform %*% diag(P[K,]) ]
  # backlog [  0  | diag(otperform)         ]
  # dimensions: n[3]*n[2]+n[1] rows(constraints), n[3]+n[1] columns(variables)
  dual.matrix=NULL
  for (assignment in c(1:n[2])){
    iden=diag(c(rep(-1,n[3])))
    performtask=ss$perform %*% diag(ss$P[assignment,])
    dual.matrix=rbind(dual.matrix,cbind(iden,performtask))
  }
  iden=matrix(c(rep(0,n[3]*n[1])),nrow=n[1])
  performtask=diag(ss$otperform)
  dual.matrix=rbind(dual.matrix,cbind(iden,performtask))

  #set up the linear programming/assignment process
  lpo=list(
    n=n,
    obj=dual.obj,
    rhs=dual.rhs,
    mat=dual.matrix,
    dir=dual.dir,
    ass=c(rep(0,n[3])),
    bl=c(rep(0,n[1])),
    istar=NULL,
    jstar=NULL,
    cost=0,
    exit=0,
    count=0)

  while (lpo$exit==0){
    lpo=iterate_dual_process(lpo)
    ss$istar[lpo$count]=lpo$istar
    ss$jstar[lpo$count]=lpo$jstar
  }
  exitmessages=c("all good","LP infeasibility","double up assignment","no reduction in demand")

  ss$solution=lpo$ass
  ss$backlog=lpo$bl
  ss$status=exitmessages[lpo$exit]
  ss
}
