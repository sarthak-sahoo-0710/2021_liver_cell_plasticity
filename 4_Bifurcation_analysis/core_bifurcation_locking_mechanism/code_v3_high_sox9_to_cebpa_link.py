import PyDSTool as dst
import numpy as np
from matplotlib import pyplot as plt

# we must give a name
DSargs = dst.args(name='Hepatoblast to Hepatocyte/Cholangiocyte cell fate decision')
# parameters
DSargs.pars = {# production rates
               'g_cebpa' : 12000, # g/k = 12000
               'g_tgfbr2': 9000, # g/k = 18000
               'g_sox9'  : 400, # g/k = 4000
               # degradation rates
               'k_cebpa' : 1.04, # hr-1
               'k_tgfbr2': 0.50, # 0.23 to 0.7 hr-1
               'k_sox9'  : 0.10, # 0.01 to 0.23 hr-1
               
               # CEBPA interactions
               'lm_cebpa_cebpa' : 2,
               'phi_cebpa_cebpa': 6000,
               'n_cebpa_cebpa'  : 3,

               'lm_cebpa_tgfbr2' : 0.1,
               'phi_cebpa_tgfbr2': 9000,
               'n_cebpa_tgfbr2'  : 3,

               'lm_cebpa_sox9' : 0.35, #base = 0.5
               'phi_cebpa_sox9': 2000,
               'n_cebpa_sox9'  : 3,
               # TGFBR2 interactions
               'lm_tgfbr2_cebpa' : 0.2,
               'phi_tgfbr2_cebpa': 6000,
               'n_tgfbr2_cebpa'  : 5,

# check
               'lm_tgfbr2_tgfbr2' : 1.25,   # for irreversibility (increasing increases the hepatocyte-hepatoblast boundary)
               'phi_tgfbr2_tgfbr2': 13000, #from 9000
               'n_tgfbr2_tgfbr2'  : 3,

# check
               'lm_tgfbr2_tgfb' : 2, #core value = 2  #from 2.25 # height of the hepatoblast-cholangiocyte sigmoid (decreaing it removes bistable region)
               'phi_tgfbr2_tgfb': 25000,
               'n_tgfbr2_tgfb'  : 3,

# check
               'lm_tgfbr2_sox9' : 0.5, #base = 0.5 # from 0.35,  # cebpa tgfbr2 bistabilty vanishes if reduced
               'phi_tgfbr2_sox9': 4000, #base = 4000# if increased to 7000 causes the upper arm to shift up by 2-3k
               'n_tgfbr2_sox9'  : 3,
               
               # SOX9 interactions
               'lm_sox9_tgfbr2' : 1.4,
               'phi_sox9_tgfbr2': 18000,
               'n_sox9_tgfbr2'  : 3,

               'lm_sox9_sox9' : 2,
               'phi_sox9_sox9': 6500,
               'n_sox9_sox9'  : 4,
               # input signal
               'tgfb' : 7600, # 7600 - 41041
               }
               
# auxiliary helper function(s) -- function name: ([func signature], definition)
DSargs.fnspecs  = {'hill': (['X','Y','lm_x_y','phi_x_y','n_x_y'], '(phi_x_y**n_x_y + lm_x_y*Y**n_x_y)/(phi_x_y**n_x_y + Y**n_x_y)') }
# rhs of the differential equation, including dummy variable w

DSargs.varspecs = {'cebpa': '(g_cebpa*\
                             hill(cebpa,cebpa,lm_cebpa_cebpa,phi_cebpa_cebpa,n_cebpa_cebpa)*\
                             hill(cebpa,tgfbr2,lm_cebpa_tgfbr2,phi_cebpa_tgfbr2,n_cebpa_tgfbr2)*\
                             hill(cebpa,sox9,lm_cebpa_sox9,phi_cebpa_sox9,n_cebpa_sox9))\
                             - k_cebpa*cebpa',
                   'tgfbr2': '(g_tgfbr2*\
                             hill(tgfbr2,tgfb,lm_tgfbr2_tgfb,phi_tgfbr2_tgfb,n_tgfbr2_tgfb)*\
                             hill(tgfbr2,tgfbr2,lm_tgfbr2_tgfbr2,phi_tgfbr2_tgfbr2,n_tgfbr2_tgfbr2)*\
                             hill(tgfbr2,sox9,lm_tgfbr2_sox9,phi_tgfbr2_sox9,n_tgfbr2_sox9)*\
                             hill(tgfbr2,cebpa,lm_tgfbr2_cebpa,phi_tgfbr2_cebpa,n_tgfbr2_cebpa))\
                             - k_tgfbr2*tgfbr2',
                   'sox9': ' (g_sox9*\
                             hill(sox9,sox9,lm_sox9_sox9,phi_sox9_sox9,n_sox9_sox9)*\
                             hill(sox9,tgfbr2,lm_sox9_tgfbr2,phi_sox9_tgfbr2,n_sox9_tgfbr2))\
                             - k_sox9*sox9'
                              }
                              
# initial conditions
DSargs.ics      = {'cebpa': 12000,'tgfbr2':2000,'sox9':5000 }

DSargs.tdomain = [0,100]                         # set the range of integration.
ode  = dst.Generator.Vode_ODEsystem(DSargs)     # an instance of the 'Generator' class.
traj = ode.compute('Model')              # integrate ODE
pts  = traj.sample(dt=0.1)                      # Data for plotting

# Prepare the system to start close to a steady state
ode.set(pars = {'tgfb': 7600} )       # Lower bound of the control parameter 'i'
ode.set(ics =  {'tgfbr2': 4000} )       # Close to one of the steady states present for i=-220

PC = dst.ContClass(ode)            # Set up continuation class

PCargs = dst.args(name='EQ1', type='EP-C')     # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
PCargs.freepars     = ['tgfb']                    # control parameter(s) (it should be among those specified in DSargs.pars)
PCargs.MaxNumPoints = 1000                      # The following 3 parameters are set after trial-and-error
PCargs.MaxStepSize  = 100
PCargs.MinStepSize  = 1
PCargs.StepSize     = 1
PCargs.LocBifPoints = 'LP'                     # detect limit points / saddle-node bifurcations
PCargs.SaveEigen    = True                     # to tell unstable from stable branches

PC.newCurve(PCargs)
PC['EQ1'].backward()
bf = list(PC['EQ1'].sol['tgfb'])
C = list(PC['EQ1'].sol['cebpa'])
T = list(PC['EQ1'].sol['tgfbr2'])
S = list(PC['EQ1'].sol['sox9'])
PC['EQ1'].forward()
bf = bf + list(PC['EQ1'].sol['tgfb'])
C = C + list(PC['EQ1'].sol['cebpa'])
T = T + list(PC['EQ1'].sol['tgfbr2'])
S = S + list(PC['EQ1'].sol['sox9'])

plt.scatter(bf,S,label="sox9",s=0.5)
plt.xlim([0,41000])
plt.legend()
plt.show()
plt.close()

plt.scatter(bf,T,label="tgfbr2",s=0.5)
plt.xlim([0,41000])
plt.legend()
plt.show()
plt.close()

plt.scatter(bf,C,label="cebpa",s=0.5)
plt.xlim([0,41000])
plt.legend()
plt.show()
plt.close()

with open("bifurcation_core_data_7_sept.txt","w") as f:
    for i,j in enumerate(bf):
        if j < 0:
            continue
        f.write(str(j)+"\t"+str(C[i])+"\t"+str(T[i])+"\t"+str(S[i])+"\n")