

import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
def readFile(path):
    f = open(path)
    first_ele = True
    for data in f.readlines():
        data = data.strip('\n')
        nums = data.split(" ")
        if first_ele:
            nums = [float(x) for x in nums ]
            matrix = np.array(nums)
            first_ele = False
        else:
            nums = [float(x) for x in nums ]
            matrix = np.c_[matrix,nums]
    return matrix
    f.close()


if __name__ == '__main__':
    
##########################################################
    
    DSCMD=np.transpose(readFile('Result_FedCET_TAU2.txt'))
    DSCMDD=np.transpose(readFile('Result_FedTrack_TAU2.txt'))
    DSCMDDD=np.transpose(readFile('Result_SCAFFOLD_TAU2.txt'))

##########################################################
    
    #DSCMD=np.transpose(readFile('Result_FedCET_TAU3.txt'))
    #DSCMDD=np.transpose(readFile('Result_FedTrack_TAU3.txt'))
    #DSCMDDD=np.transpose(readFile('Result_SCAFFOLD_TAU3.txt'))
    
    

##########################################################
    
    #DSCMD=np.transpose(readFile('Result_FedCET_TAU4.txt'))
    #DSCMDD=np.transpose(readFile('Result_FedTrack_TAU4.txt'))
    #DSCMDDD=np.transpose(readFile('Result_SCAFFOLD_TAU4.txt'))


##########################################################
    
    #DSCMD=np.transpose(readFile('Result_FedCET_TAU5.txt'))
    #DSCMDD=np.transpose(readFile('Result_FedTrack_TAU5.txt'))
    #DSCMDDD=np.transpose(readFile('Result_SCAFFOLD_TAU5.txt'))    
    
    
    DSCMD1=DSCMD
    DSCMDD1=DSCMDD
    DSCMDDD1=DSCMDDD
    
    
    
    
    k=440
    kk=61
    
    kk_OTHER_FedCET=31
    
    BASE1=np.zeros(k)
    BASE2=np.zeros(k)
    BASE3=np.zeros(kk)
    BASE4=np.zeros(kk)
    BASE5=np.zeros(k)
    BASE6=np.zeros(k)
    BASE7=np.zeros(kk_OTHER_FedCET)
    BASE8=np.zeros(kk_OTHER_FedCET)
    BASE9=np.zeros(k)
    BASE10=np.zeros(k)
    BASE11=np.zeros(kk_OTHER_FedCET)
    BASE12=np.zeros(kk_OTHER_FedCET)
    
    
    
    
    k0=0
    
    
    for j in range(0,k-1):
        BASE1[j]=j
        BASE2[j]=DSCMD1[j+1]/math.pow(10,14)
        
        
    for j in range(0,k-1):
        BASE5[j]=j
        BASE6[j]=DSCMDD1[j+1]/math.pow(10,14)
        
    for j in range(0,k-1):
        BASE9[j]=j
        BASE10[j]=DSCMDDD1[j+1]/math.pow(10,14)
    
    
    
    for jj in range(0,kk):
        if jj%2==0:
            add=jj#int(jj/2)
        else:
            add=jj #int((jj+1)/2)
        BASE3[jj]=BASE1[k0+add]
        BASE4[jj]=BASE2[k0+jj]
      
        
    for jj in range(0,kk_OTHER_FedCET):
        BASE7[jj]=BASE5[k0+2*jj]
        BASE8[jj]=BASE6[k0+jj]
       
    
    
    for jj in range(0,kk_OTHER_FedCET):
        BASE11[jj]=BASE9[k0+2*jj]
        BASE12[jj]=BASE10[k0+jj]
        
        
    
    x = BASE3
    y = BASE4
    
    
    x2=BASE7
    y2=BASE8
    
    x3=BASE11
    y3=BASE12
    
    
   
    
    
    
    
    plt.xticks(fontsize = 10, fontname = 'times new roman')
    plt.yticks(fontsize = 10, fontname = 'times new roman')
 
    
    plt.semilogy(x3,y3,label='SCAFFOLD',color='b',linewidth=2.4,linestyle='--')
    plt.semilogy(x2,y2,label='FedTrack',color='r',linewidth=2.4,linestyle='-.')
    plt.semilogy(x,y,label='FedCET', color='y',linewidth=2.4,linestyle='solid')
    
    
    font1 = {'family' : 'Times New Roman','weight' : 'normal','size'   : 18, }
 
    plt.legend(prop = {'size':13}) 
    plt.grid()
    
    plt.xlabel('Number of Shared Messages',fontsize=13)
    plt.ylabel('Convergence Errors',fontsize=13)
    
    
    
    plt.savefig(fname="different_algorithm_TAU2.jpg",format="jpg", bbox_inches='tight')
    plt.savefig(fname="different_algorithm_TAU2.pdf",format="pdf", bbox_inches='tight')
    plt.show() 
    
    
 #   plt.savefig(fname="different_algorithm_TAU3.jpg",format="jpg", bbox_inches='tight')
 #   plt.savefig(fname="different_algorithm_TAU3.pdf",format="pdf", bbox_inches='tight')
 #   plt.show() 
    
    
 #   plt.savefig(fname="different_algorithm_TAU4.jpg",format="jpg", bbox_inches='tight')
 #   plt.savefig(fname="different_algorithm_TAU4.pdf",format="pdf", bbox_inches='tight')
 #   plt.show() 
    
    
 #   plt.savefig(fname="different_algorithm_TAU5.jpg",format="jpg", bbox_inches='tight')
 #   plt.savefig(fname="different_algorithm_TAU5.pdf",format="pdf", bbox_inches='tight')
 #   plt.show() 
    
 
