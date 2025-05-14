

import numpy as np
import random 
import matplotlib.pyplot as plt
import math
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
   

#tau=2
#tau=3
#tau=4
tau=5
  
iter=550
iter_update=iter*tau
N=10
N_row=60
N_col=60
n=60
ni=10

b=np.zeros((N,ni,n))
C=np.zeros((N,n,n))
D=np.zeros((N,n))

for s1 in range(0,N):
    for s2 in range(0,ni):
        for s3 in range(0,n):
            b[s1][s2][s3]=random.uniform(-10,10)


         
    
mu=2  
L_Max=1.5

number_count=1

step_size_min=mu*mu*(1/4)
if step_size_min>1/(10*tau*tau*math.pow(L_Max,4)):
    step_size_min=1/(math.pow(1+2/tau,tau+2)*tau*tau*math.pow(L_Max,4))

step_size_max=2/(tau*L_Max)

distance_each=(step_size_max-step_size_min)/100

step_size=step_size_min


h=0.001*step_size_min

xi=tau*mu

K=math.ceil((step_size_max-step_size_min)/h)

for k in range(0,K):
    print(k)
    parameter1=2*tau*mu*step_size-xi*step_size+math.pow(tau,4)*math.pow(1+2/tau,2*tau-2)*math.pow(step_size,4)*math.pow(L_Max,4)-math.pow(tau,4)*math.pow(1+2/tau,2*tau-2)*math.pow(step_size,3)*math.pow(L_Max,4)*(2/(xi))-math.pow(tau,2)*L_Max*mu*math.pow(step_size,2)
    parameter2=1-xi*step_size-((2/(xi*step_size))-1)*math.pow(tau,2)*math.pow(1+2/tau,2*tau-2)*math.pow(step_size,2)*math.pow(L_Max,2)
    if parameter1<0:
        step_size=step_size
    else:
        if parameter2<0:
            step_size=step_size
        else:
            step_size=step_size+h
            number_count=1+number_count
        


c=mu/(2*mu*step_size+8)

    

C=np.zeros(n)

B=np.zeros((N,n))

for s6 in range(0,N):
    for s7 in range(0,ni):
        B[s6]=B[s6]+b[s6][s7]
    B[s6]=B[s6]/ni  
        
        

for s4 in range(0,N):
    for s5 in range(0,ni):
        C=C+b[s4][s5]
    
X_OPT=(1/(2*ni*N))*C



####################################################################


def obj(XX):
    result=0
    for p in range(0,N):
        for q in range(0,ni):
            result=result+np.linalg.norm(XX-b[p][q])*np.linalg.norm(XX-b[p][q])
    final_result=result/(2*N*ni)+np.linalg.norm(XX)*np.linalg.norm(XX)
    return final_result

####################################################################

Result_FedCET=np.zeros(iter)

number_FedCET=0



X_my1=0*np.zeros((iter_update,N,n))


for j in range(0,N):
    X_my1[1,j]=X_my1[0,j]-step_size*(4*X_my1[0,j]-2*B[j])
    
for i in range(0,N):
    VALUE=0*X_my1[0,1]
    for j in range(0,N):
        g1=4*X_my1[1,j]-2*B[j] 
        g2=4*X_my1[0,j]-2*B[j] 
        VALUE=VALUE+(2*X_my1[1,j]-X_my1[0,j]-step_size*g1+step_size*g2)
    VALUE=VALUE/N 
    X_my1[2,i]=c*step_size*VALUE+(1-c*step_size)*(2*X_my1[1,i]-X_my1[0,i]-step_size*(4*X_my1[1,i]-2*B[i])+step_size*(4*X_my1[0,i]-2*B[i]))
    
value_result1=X_my1[2,0]
for i in range(1,N):
    value_result1=value_result1+X_my1[2,i]
value_result1=value_result1/N
print(np.linalg.norm(value_result1-X_OPT))

Result_FedCET[number_FedCET]=np.linalg.norm(value_result1-X_OPT)
number_FedCET=number_FedCET+1
    
for t in range(2,iter_update-5):  
    for i in range(0,N):
        if t%tau==0:
            VALUE=0*X_my1[t,i]
            for j in range(0,N):
                g1=4*X_my1[t,j]-2*B[j]  
                g2=4*X_my1[t-1,j]-2*B[j] 
                VALUE=VALUE+(2*X_my1[t,j]-X_my1[t-1,j]-step_size*g1+step_size*g2)
            VALUE=VALUE/N 
            X_my1[t+1,i]=c*step_size*VALUE+(1-c*step_size)*(2*X_my1[t,i]-X_my1[t-1,i]-step_size*(4*X_my1[t,i]-2*B[i])+step_size*(4*X_my1[t-1,i]-2*B[i]))
        else:
            g1=4*X_my1[t,i]-2*B[i]  
            g2=4*X_my1[t-1,i]-2*B[i] 
            X_my1[t+1,i]=2*X_my1[t,i]-X_my1[t-1,i]-step_size*g1+step_size*g2
    
    if t%tau==0:
        print(t)
        value_result1=X_my1[t,0]
        for i in range(1,N):
            value_result1=value_result1+X_my1[t,i]
        value_result1=value_result1/N
        print(np.linalg.norm(value_result1-X_OPT))
        Result_FedCET[number_FedCET]=np.linalg.norm(value_result1-X_OPT)*math.pow(10,14)
        number_FedCET=number_FedCET+1
        
        
        
##########################################################
Result_SCAFFOLD=np.zeros(iter)

X_my2=0*np.zeros((N,n))
Y2=0*np.zeros((N,n))  
X2=0*np.zeros(n)  
step_size2g=1
step_size2l=1/(81*L_Max*tau)
C2=0*np.zeros((N,n))
C_PLUS=0*np.zeros((N,n))
C_ALL=0*np.zeros(n)  

Delta_Y2=0*np.zeros((N,n))  
Delta_C2=0*np.zeros((N,n))

for t in range(0,iter):
    for i in range(0,N):
        Y2[i]=X2
        for k in range(0,tau):
            Y2[i]=Y2[i]-step_size2l*((4*Y2[i]-2*B[i])-C2[i]+C_ALL)
        C_PLUS[i]=4*X2-2*B[i]
        Delta_Y2[i]=Y2[i]-X2
        Delta_C2[i]=C_PLUS[i]-C2[i]
        C2[i]=C_PLUS[i]
    AAA=0*Delta_Y2[0]
    BBB=0*Delta_C2[0]
    for j in range(0,N):
        AAA=AAA+(1/N)*Delta_Y2[j]
        BBB=BBB+(1/N)*Delta_C2[j]
    X2=X2+step_size2g*AAA
    C_ALL=C_ALL+BBB
    print(t)
    print(np.linalg.norm(X2-X_OPT))
    Result_SCAFFOLD[t]=np.linalg.norm(X2-X_OPT)*math.pow(10,14)
     
            


##########################################################
Result_FedTrack=np.zeros(iter)

step_size3=1/(18*L_Max*tau)


X_my3=0*np.zeros((N,n))


GG_my3=np.zeros(n)
GG_AVG_my3=np.zeros(n)
ave_X_my3=np.zeros(n)

for t in range(0,iter):
    
    
    
    
    for i in range(0,N):  
        X_my3[i]=ave_X_my3
        GG_FITST_my3=4*X_my3[i]-2*B[i]   # np.dot(np.dot(np.transpose(A[i]),A[i]),X_my3[i])-np.dot(np.transpose(A[i]),b[i])
        
        for k in range(0,tau):
            X_my3[i]=X_my3[i]-step_size3*(GG_AVG_my3-GG_FITST_my3+(4*X_my3[i]-2*B[i]))
            
            
   
    ave_X_my3=np.zeros(n)
    GG_my3=np.zeros(n)
    
    for j in range(0,N):
        ave_X_my3=ave_X_my3+X_my3[j]
        
    ave_X_my3=ave_X_my3/N
    
    for ii in range(0,N): 
        GG_my3=GG_my3+4*ave_X_my3-2*B[ii]  # np.dot(np.dot(np.transpose(A[ii]),A[ii]),ave_X_my3)-np.dot(np.transpose(A[ii]),b[ii])
        
    
    GG_AVG_my3=GG_my3/N
    
    
    print(t)
    
    value_result3=X_my3[0]
    for i in range(1,N):
        value_result3=value_result3+X_my3[i]
    value_result3=value_result3/N
    
    print(np.linalg.norm(value_result3-X_OPT))
    Result_FedTrack[t]=np.linalg.norm(value_result3-X_OPT)*math.pow(10,14)
    
################################################################


#np.savetxt('Result_FedCET_TAU2.txt',Result_FedCET,fmt="%f",delimiter=" ") 
#np.savetxt('Result_FedTrack_TAU2.txt',Result_FedTrack,fmt="%f",delimiter=" ")
#np.savetxt('Result_SCAFFOLD_TAU2.txt',Result_SCAFFOLD,fmt="%f",delimiter=" ")


################################################################


#np.savetxt('Result_FedCET_TAU3.txt',Result_FedCET,fmt="%f",delimiter=" ") 
#np.savetxt('Result_FedTrack_TAU3.txt',Result_FedTrack,fmt="%f",delimiter=" ")
#np.savetxt('Result_SCAFFOLD_TAU3.txt',Result_SCAFFOLD,fmt="%f",delimiter=" ")



################################################################


#np.savetxt('Result_FedCET_TAU4.txt',Result_FedCET,fmt="%f",delimiter=" ") 
#np.savetxt('Result_FedTrack_TAU4.txt',Result_FedTrack,fmt="%f",delimiter=" ")
#np.savetxt('Result_SCAFFOLD_TAU4.txt',Result_SCAFFOLD,fmt="%f",delimiter=" ")



################################################################


np.savetxt('Result_FedCET_TAU5.txt',Result_FedCET,fmt="%f",delimiter=" ") 
np.savetxt('Result_FedTrack_TAU5.txt',Result_FedTrack,fmt="%f",delimiter=" ")
np.savetxt('Result_SCAFFOLD_TAU5.txt',Result_SCAFFOLD,fmt="%f",delimiter=" ")