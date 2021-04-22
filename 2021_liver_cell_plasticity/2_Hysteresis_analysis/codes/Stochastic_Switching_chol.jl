using DifferentialEquations
using Plots
pyplot()

function hill(X,Y,lm_x_y,phi_x_y,n_x_y)::Float64
    h = (big(phi_x_y)^n_x_y + big(lm_x_y)*big(Y)^n_x_y)/(big(phi_x_y)^n_x_y + big(Y)^n_x_y)
    return h
end

function get_parameters()
    #production rates
    g_cebpa = 12000
    g_tgfbr2 = 9000
    g_sox9 = 400

    #degradation rates
    k_cebpa = 1.04
    k_tgfbr2 = 0.50
    k_sox9 = 0.10

    #CEBPA interactions
    lm_cebpa_cebpa = 2
    phi_cebpa_cebpa = 6000
    n_cebpa_cebpa = 3

    lm_cebpa_tgfbr2 = 0.1
    phi_cebpa_tgfbr2 = 9000
    n_cebpa_tgfbr2 = 3

    lm_cebpa_sox9 = 0.5
    phi_cebpa_sox9 = 2000
    n_cebpa_sox9 = 3

    #TGFBR2 interactions
    lm_tgfbr2_cebpa = 0.2
    phi_tgfbr2_cebpa= 6000
    n_tgfbr2_cebpa = 5

    lm_tgfbr2_tgfbr2 = 1.25   # for irreversibility (increasing increases the hepatocyte-hepatoblast boundary)
    phi_tgfbr2_tgfbr2= 13000
    n_tgfbr2_tgfbr2  = 3

    lm_tgfbr2_tgfb = 2   # height of the hepatoblast-cholangiocyte sigmoid (decreaing it removes bistable region)
    phi_tgfbr2_tgfb= 25000
    n_tgfbr2_tgfb  = 3

    lm_tgfbr2_sox9 = 0.5  # cebpa tgfbr2 bistabilty vanishes if reduced
    phi_tgfbr2_sox9= 4000
    n_tgfbr2_sox9  = 3
    # SOX9 interactions
    lm_sox9_tgfbr2 = 1.4
    phi_sox9_tgfbr2= 18000
    n_sox9_tgfbr2  = 3

    lm_sox9_sox9 = 2
    phi_sox9_sox9= 6500
    n_sox9_sox9  = 4
    # input signal
    tgfb = 20000 # 7600 - 41041
    
    return g_cebpa, g_tgfbr2, g_sox9, k_cebpa, k_tgfbr2, k_sox9, lm_cebpa_cebpa, phi_cebpa_cebpa,
        n_cebpa_cebpa, lm_cebpa_tgfbr2, phi_cebpa_tgfbr2, n_cebpa_tgfbr2,
        lm_cebpa_sox9, phi_cebpa_sox9, n_cebpa_sox9, lm_tgfbr2_cebpa, phi_tgfbr2_cebpa,
        n_tgfbr2_cebpa, lm_tgfbr2_tgfbr2, phi_tgfbr2_tgfbr2, n_tgfbr2_tgfbr2,
        lm_tgfbr2_tgfb, phi_tgfbr2_tgfb, n_tgfbr2_tgfb, lm_tgfbr2_sox9, phi_tgfbr2_sox9,
        n_tgfbr2_sox9, lm_sox9_tgfbr2, phi_sox9_tgfbr2, n_sox9_tgfbr2,
        lm_sox9_sox9, phi_sox9_sox9, n_sox9_sox9, tgfb
    
end

function deterministic(du,u,p,t)
    g_cebpa, g_tgfbr2, g_sox9, k_cebpa, k_tgfbr2, k_sox9, lm_cebpa_cebpa, phi_cebpa_cebpa,
        n_cebpa_cebpa, lm_cebpa_tgfbr2, phi_cebpa_tgfbr2, n_cebpa_tgfbr2,
        lm_cebpa_sox9, phi_cebpa_sox9, n_cebpa_sox9, lm_tgfbr2_cebpa, phi_tgfbr2_cebpa,
        n_tgfbr2_cebpa, lm_tgfbr2_tgfbr2, phi_tgfbr2_tgfbr2, n_tgfbr2_tgfbr2,
        lm_tgfbr2_tgfb, phi_tgfbr2_tgfb, n_tgfbr2_tgfb, lm_tgfbr2_sox9, phi_tgfbr2_sox9,
        n_tgfbr2_sox9, lm_sox9_tgfbr2, phi_sox9_tgfbr2, n_sox9_tgfbr2,
        lm_sox9_sox9, phi_sox9_sox9, n_sox9_sox9, tgfb = get_parameters()
    du[1] = (g_cebpa * 
            hill(u[1],u[1],lm_cebpa_cebpa,phi_cebpa_cebpa,n_cebpa_cebpa) * 
            hill(u[1],u[2],lm_cebpa_tgfbr2,phi_cebpa_tgfbr2,n_cebpa_tgfbr2)*
            hill(u[1],u[3],lm_cebpa_sox9,phi_cebpa_sox9,n_cebpa_sox9)) - 
            k_cebpa*u[1]
    du[2] = (g_tgfbr2*
            hill(u[2],tgfb,lm_tgfbr2_tgfb,phi_tgfbr2_tgfb,n_tgfbr2_tgfb)*
            hill(u[2],u[2],lm_tgfbr2_tgfbr2,phi_tgfbr2_tgfbr2,n_tgfbr2_tgfbr2)*
            hill(u[2],u[3],lm_tgfbr2_sox9,phi_tgfbr2_sox9,n_tgfbr2_sox9)*
            hill(u[2],u[1],lm_tgfbr2_cebpa,phi_tgfbr2_cebpa,n_tgfbr2_cebpa)) -
            k_tgfbr2*u[2]
    du[3] = (g_sox9*
            hill(u[3],u[3],lm_sox9_sox9,phi_sox9_sox9,n_sox9_sox9)*
            hill(u[3],u[2],lm_sox9_tgfbr2,phi_sox9_tgfbr2,n_sox9_tgfbr2)) - 
            k_sox9*u[3]
end 

function noise(du,u,p,t)
    
    du[1] = 1100
    du[2] = 1100
    du[3] = 110
    
end







function noise(du,u,p,t)
    
    du[1] = 1100
    du[2] = 1100
    du[3] = 110
    
end

u0_cebpa = 11000#2100
u0_tgfbr2 = 2500#12500
u0_sox9 = 5100#6500
u0 = [u0_cebpa,u0_tgfbr2,u0_sox9]

tspan = (0.0,100.0)

f = open("backward_20k.txt","w")
for x in 1:300
    sde = SDEProblem(deterministic,noise,u0,tspan)
    sol = solve(sde);
    #println(last(sol)[1])
    write(f,string(last(sol)[1]),"\t",string(last(sol)[2]),"\t",string(last(sol)[3]),"\n")
end
close(f)




