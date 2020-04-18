import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt







def sigmatot(ymin,ymax,m,R1,s12,alpha_s,RG,s2):
    #Total cross section
    I=quad(dsigmadyaux,ymin,ymax,args=(m,R1,s12,alpha_s,RG,s2))
    return I[0]

def dsigmadyaux(y,m,R1,s12,alpha_s,RG,s2):
    #Auxiliary function to be used in sigmatot, that is, to be integrated over y (rapidity).
    res=dsigmady(m,R1,y,s12,alpha_s,RG,s2)
    return res
        

def dsigmady(m,R1,y,s12,alpha_s,RG,s2):    
    #Cross section distribution in rapidity 
    s=np.power(s12,2.0)
    b=4.0
    fac=1e6*0.389*s2
    ct=(fac*np.power(np.pi,4.0)*m*R1*np.power(alpha_s,2.0))/(b*b)
    I=quad(func_auxdsdy,0.0, np.inf,args=(m,y,alpha_s,s))
    res=ct*np.power(I[0],2.0)*np.power(RG,4.0)
    return res
        
def func_auxdsdy(qt,m,y,alpha_s,s):
    #Auxiliary function to be used in dsigmady, that is, to be integrated over Q.
    qt2=qt*qt
    x1=(m/np.sqrt(s))*np.exp(y)
    x2=(m/np.sqrt(s))*np.exp(-y)
    a1=np.power(1.0/qt2,2.0)
    b1=ugd_gbw(x1,qt2,alpha_s)*ugd_gbw(x2,qt2,alpha_s)
    c1=(2.0*qt)/((m**2.0 + 1.0*qt2)**2.0)
    res=a1*b1*c1
    return res


def ugd_gbw(x,kt2,alpha_s):
    #unintegrated gluon distribution based on GBW dipole model.
    sigma0=27.32/0.389
    x0=4.2e-5
    lamb=0.248
    ro2=np.power(x/x0,lamb)
    fac=np.power((1.0-x),5.0)
    res=(3.0*sigma0*ro2*np.power(kt2,2.0)*np.exp(-ro2*kt2)*fac)/(4.0*np.power(np.pi,2.0)*alpha_s)
    return res



mass=3.414
rp2=0.075/1.45
w=1960.0
s2=0.058
ymin=-5.0
ymax=5.0
alphas=0.335   




def totalcross():
    print('---------------------Total Cross Section-----------------------------')
    print("Total cross section ", sigmatot(ymin,ymax,mass,rp2,w,alphas,1.0,s2) ,"nb" )
    print('---------------------------------------------------------------------')
    print("")



xdata=[]
ydata=[]
def rapidity():
    #print('---------------------Dsigma/dy [nb]----------------------------------')
    #print('----Y ---Dsigma/dy [nb]----------------------------------------------')
    for y in np.arange(0.0,5.0,0.1):
        res=dsigmady(mass,rp2,y,w,alphas,1.0,s2)
        xdata.append(y)
        ydata.append(res)
        #print( y , "  " ,  res)  
        #print('---------------------------------------------------------------------')


xugd=[]
yugd=[]
def ugd():
    for kt2 in np.arange(0.0,15.0,0.1):
        x=mass/14000.0
        res=ugd_gbw(x,kt2,alphas)*alphas
        xugd.append(kt2)
        yugd.append(res)
        
    



def plots(op):
     if(op==1):
         plt.figure(1)
         #plt.xscale('log')
         plt.yscale('log')
         plt.plot(xdata,ydata,'g',label=r'$GBW$',linewidth=1.5)
         plt.title('',fontsize=15)
         plt.xlabel(r'$Y$',fontsize=15)
         plt.ylabel(r'$d\sigma/dy  \ [nb]$',fontsize=15)
         plt.legend(loc='upper right')
         plt.xticks(np.arange(0.0,5.5,0.5))
        #plt.yticks(np.arange(1,2.2,0.1))
        #plt.text(0.02,1.3, r'an equation: $E=mc^2$', fontsize=15)
         #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #plt.xlim(34,10000)
         plt.ylim(1,1e3)
         #plt.savefig('pp-yns-xg-w-0.4.pdf')  
         plt.show()
     if(op==2):
         plt.figure(1)
         #plt.xscale('log')
         #plt.yscale('log')
         plt.plot(xugd,yugd,'g',label=r'$GBW$',linewidth=1.5)
         plt.title('',fontsize=15)
         plt.xlabel(r'$Q_{T}^{2} \ GeV^{2}$ ',fontsize=15)
         plt.ylabel(r'$\alpha_{s}(m^{2}) \ fg(x=m/\sqrt{s},Q_{T}^{2},m^{2})$ ',fontsize=15)
         plt.legend(loc='upper right')
         plt.xticks(np.arange(0.0,15,2.5))
         plt.yticks(np.arange(0.0,4.0,0.5))
        #plt.text(0.02,1.3, r'an equation: $E=mc^2$', fontsize=15)
         #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #plt.xlim(34,10000)
         #plt.ylim(1,1e3)
         #plt.savefig('pp-yns-xg-w-0.4.pdf')  
         plt.show()

#rapidity()
ugd()
plots(2)

