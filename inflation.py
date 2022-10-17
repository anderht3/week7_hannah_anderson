#initial equations that are needed: 
#7.14 - phi** + 3Hphi* + v = 0
#7.8 - rho = 1/2 phi*^2 + v
#Friedman - H^2 = 8piG rho/3


from scipy.integrate import solve_ivp

def first_deg(t, y): return(y)

sol = solve_ivp(second_deg, [0,50],[5])
print(sol.t)
print(sol.y)

#def second_deg(t,y): return
