library(matlib)


#probability matrix for first action
p_first = matrix(c(0),91,91)
p_first[,1]=seq(0.1,1,0.01)

for(i in (1:90)){
  p_first[i,1+i]=1-p_first[i,1]
}
  
#stationary distribution
pi = matrix(c(0),1,91)
pi[1]=1
convergence = 1
while(convergence>0.000001){
  pi_new = pi %*% p_first
  pi_dif = pi_new - pi
  convergence = max(pi_dif)
  pi = pi_new
}

r_first = matrix(c(seq(0.1,1,0.01)),91,1)
avg_replacement = pi %*% r_first

#Poisson function
Poisson <- function(r,p){
  coefs = matrix(c(0),92,92)
  coefs[1:91,92] = 1
  coefs[92,1]=1
  
  coef_result = matrix(c(0),92,1)
  
  for (i in (1:91)){
    coefs[i,1:91] =  - p[i,1:91]
    coefs[i,i] = 1 - p[i,i]
    coef_result[i,1] = r[i,1]
  }
  
  A = solve(coefs,coef_result,fractions=TRUE)
  
  ret = list("V"=A[1:91],"g"=A[92])
  return(ret)
}


##poisson for first action
poisson_first = Poisson(r_first,p_first)
g_first_poisson = poisson_first$g
V_first_poisson = poisson_first$V

#preventive Replacement cost and probability matrix
p_second = matrix(c(0),91,91)
p_second[1:91,1] = 1

r_second = matrix(c(0.5),91,1)


###policy iteration
##initial policy is no preventive replacement at all (1st action for each state)
R_initial = matrix(c(1),91,1)

R = R_initial

change = 1000
while(change!=0){
  r_policy_iteration = matrix(c(0),91,1)
  p_policy_iteration = matrix(c(0),91,91)
  for(i in (1:91)){
    r_policy_iteration[i,1] = ifelse(R[i,1]==1,r_first[i,1],r_second[i,1])
    p_policy_iteration[i,1:91] = if(R[i,1]==1) p_first[i,1:91] else p_second[i,1:91]
  }
  poisson_policy_iteration = Poisson(r_policy_iteration,p_policy_iteration)
  first_action = (r_first + p_first %*% poisson_policy_iteration$V)
  second_action = (r_second + p_second %*% poisson_policy_iteration$V)
  R_new = ifelse(first_action >second_action ,2,1)
  change = sum(abs(R_new - R))
  print(change)
  R  = R_new
}

g_policy_iteration = (pmin(first_action,second_action) - poisson_policy_iteration$V)[1]
R_optimal_policy = R


##value iteration
##initial value matrix is all zeros
V_initial = matrix(c(0),91,1)

V = V_initial
change = 100
while(change!=1){
  first_action = (r_first + p_first %*% V)
  second_action = (r_second + p_second %*% V)
  V_new = pmin(first_action,second_action)
  V_dif = V_new - V
  change = length(unique(round(V_dif,6)))
  V = V_new
}
g_value_iteration = V_dif[1]
R_optimal_value = ifelse(first_action >second_action ,2,1)


