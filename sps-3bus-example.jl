using JuMP, GLPK,  LinearAlgebra, DataFrames
model = Model(GLPK.Optimizer)

#Data
L = [150,200,300]
P = [175,450,25]

cdP = [10,10,10]
cdL = [20,20,20]

dLmax = [0,0,0]
dLmin = -L

dPmax = [200,200,200]
dPmin = -P

#reactances in contingency
X12 = 0.05
X13 = 0.05
X23 = 0.1
#max capacity of Transmission lines
F12max = 200
F13max = 200
F23max = 100

#Variables
@variable(model,dP[1:3])
@variable(model,dL[1:3])
@variable(model,theta[1:3])

#flows
@variable(model,F12)
@variable(model,F13)
@variable(model,F23)

#binary
@variable(model,bdL[1:3],Bin)
@variable(model,bdP[1:3],Bin)

#Transmission constraints
@constraint(model,cons_fmax12a,F12<=F12max)
@constraint(model,cons_fmax12b,F12>=-F12max)

@constraint(model,cons_fmax13a,F13<=F13max)
@constraint(model,cons_fmax13b,F13>=-F13max)

@constraint(model,cons_fmax23a,F23<=F23max)
@constraint(model,cons_fmax23b,F23>=-F23max)

#Angles
@constraint(model,F12==(theta[2]-theta[1])/X12)
@constraint(model,F23==(theta[3]-theta[2])/X23)
@constraint(model,F13==(theta[3]-theta[1])/X13)


#Nodal balacing in contingency
@constraint(model,cons_balance_1,P[1]+dP[1]==L[1]+dL[1]+F12+F13)
@constraint(model,cons_balance_2,P[2]+dP[2]+F12==L[2]+dL[2]+F23)
@constraint(model,cons_balance_3,P[3]+dP[3]+F13+F23==L[3]+dL[3])

#limits of corrective actions
@constraint(model,[i=1:3],dL[i]>=dLmin[i]*bdL[i])
@constraint(model,[i=1:3],dL[i]<=dLmax[i]*bdL[i])

@constraint(model,[i=1:3],dP[i]>=dPmin[i]*bdP[i])
@constraint(model,[i=1:3],dP[i]<=dPmax[i]*bdP[i])


#set objective function
@objective(model,Min,dot(bdP,cdP)+dot(bdL,cdL))

JuMP.optimize!(model)
println(termination_status(model))

println("Powers")
for i in 1:3
    println("P[",i,"] = ",P[i], " \t dP[",i,"] = ",value(dP[i]))
end

println("Loads")
for i in 1:3
    println("L[",i,"] = ",L[i], " \t dL[",i,"] = ",value(dL[i]))
end

println("angles")
for i in 1:3
    println("theta[",i,"] = ",value(theta[i]))
end

println("Flows")
println("F12 = ",value(F12))
println("F13 = ",value(F13))
println("F123 = ",value(F23))

println("Objective = ",objective_value(model))
