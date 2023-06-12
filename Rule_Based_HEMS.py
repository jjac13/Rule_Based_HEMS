# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################
import time
from numpy import append, array, cos, pi, arange

###########################   Functions definition  ###########################

def heat_loss(U, Area, T_in, T_out):
    return U * Area * (T_out - T_in)
    
def new_house_Temperature(T_0, Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_HP = 0, Qdot_TESS = 0, Qdot_SC = 0, dt=1):
    T_new = T_0 + dt*(Qdot_roof + Qdot_windows + Qdot_wall + Qdot_HP + Qdot_TESS + Qdot_SC)/(mc_T)
    
    return T_new

def update_TESS(active, T_0, T_soil, Qdot_SC = 0, Qdot_HP = 0, Tmax = 95 + 273, Tmin = 50 + 273, mdot = 0.1, T_network = 40 + 273, Qdot_SD = 100, efficiency = 0.8, m = 4000, c = 4200, dt=1*3600):

    if active and T_0 >= Tmin:
        Qdot_TESS = mdot*c*(T_0*0.97 - T_network)
    else:
        Qdot_TESS = 0

    
    if T_0 <= Tmin:         # TESS discharged, only charge
        if T_0 <= T_soil:   # Code for heat dissipation/absoprtion through the soil needed
            T_new = T_0
            
        else:
            T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
            
    elif T_0 <= Tmax:       # TESS available for charge and discharge
        
        T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
        
    else:                   # TESS fully charged, only discharge
        
        T_new = T_0 + (- Qdot_SD - Qdot_TESS)*dt/(m*c)        
        
    return [T_new, Qdot_TESS*efficiency]          

def Qdot_SolarCollector(active, G, A = 6, SC_eff = 0.45, dt = 1*3600):
    
    if active:
        return A*SC_eff*G/dt
    else:
        return 0


def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
    
def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)

# Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet
def update_BESS(SoC_0, P_Load, P_PV, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):
    
    
    E_BESS_0 = Capacity_BESS*SoC_0
    
    if SoC_0 >= SoCmax:                  # Battery can only discharge

        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
            
#            state = 1

        
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max  <= BESS_perm_max(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = P_Grid_max
                
#                state = 2
                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0) - P_PV + P_Load     
                
#                state = 3
    
    elif SoC_0 > SoCmin:                  # Battery can charge and discharge
        
        if (P_Load-P_PV) <= P_Grid_max:        # No peakshaving needed
            
            if P_Load >= P_PV:                    # PV below demanded load
                P_BESS = 0
                E_BESS = E_BESS_0*(1-P_SD)
                SoC_BESS = E_BESS/Capacity_BESS
                P_Grid = P_Load - P_PV     

#                state = 4

            else:                                                   # Surplus of PV power
                if P_PV - P_Load <= -BESS_perm_min(SoC_0):    # Below the BESS max power
                    P_BESS = P_Load - P_PV
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS /  Capacity_BESS
                    P_Grid = 0
                    
#                    state = 5

                    
                else:                                                       # Above the BESS max power
                    P_BESS = BESS_perm_min(SoC_0)
                    E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                    E_BESS = E_BESS*(1-P_SD)
                    SoC_BESS = E_BESS /  Capacity_BESS
                    P_Grid = -BESS_perm_min(SoC_0) - P_PV + P_Load
                    
#                    state = 6

                    
        else:                                                       # Peakshaving needed
            if P_Load - P_PV - P_Grid_max <= BESS_perm_max(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV - P_Grid_max
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = P_Grid_max
                
#                state = 7

                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_max(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_max(SoC_0) - P_PV + P_Load  
                
#                state = 8

    
    else: # self.BESS.SoC[1,i] <= SoCmin:                 # Battery can only charge     
        
        if P_Load >= P_PV:                        # PV below demanded load
            P_BESS = 0
            E_BESS = E_BESS_0*(1-P_SD)
            SoC_BESS = E_BESS/Capacity_BESS
            P_Grid = P_Load - P_PV 
            
#            state = 9
             
        else:                                                       # Surplus of PV power
            
            if P_PV - P_Load <= -BESS_perm_min(SoC_0):    # Below the BESS max power
                P_BESS = P_Load - P_PV
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = 0
                
#                state = 10

                
                
            else:                                                       # Above the BESS max power
                P_BESS = BESS_perm_min(SoC_0)
                E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
                E_BESS = E_BESS*(1-P_SD)
                SoC_BESS = E_BESS /  Capacity_BESS
                P_Grid = -BESS_perm_min(SoC_0) - P_PV + P_Load      
                
#                state = 11
                
                
    return [SoC_BESS, P_Grid, P_BESS]

def HP_Power(active, P_in = 2.7, COP = 4.1):
    if active:
        return [P_in, P_in*COP*1000]
    else:
        return [0,0]

###############################################################################
##################################   Example ##################################

# Convective heat transfer coefficients [W/m^2K]
h_air_wall   = 0.9#24       # Indoor air -> walls, scaled to R-value of a C-label house
h_wall_atm   = 0.9#34       # Walls -> atmosphere, scaled to R-value of a C-label house
h_air_window = 25            # Indoor air -> windows
h_window_atm = 32            # Windows -> atmosphere
h_air_roof   = 12            # Indoor air -> roof
h_roof_atm   = 38            # Roof -> atmosphere

## House

# Air
c_air        = 1005.4        # Specific heat of air at 273 K [J/kgK]
airDensity   = 1.025         # Densiity of air at 293 K [kg/m^3]
kAir         = 0.0257        # Thermal conductivity of air at 293 K [W/mK]


# Windows (glass)
n1_window     = 3            # Number of windows in room 1
n2_window     = 2            # Number of windows in room 2
n3_window     = 2            # Number of windows in room 3
n4_window     = 1            # Number of windows in room 4
htWindows     = 1            # Height of windows [m]
widWindows    = 1            # Width of windows [m]
windows_area  = (n1_window + n2_window + n3_window + n4_window) * htWindows * widWindows
LWindow       = 0.004        # Thickness of a single window pane [m]
LCavity       = 0.014        # Thickness of the cavity between the double glass window [m]  
windowDensity = 2500         # Density of glass [kg/m^3]
c_window      = 840          # Specific heat of glass [J/kgK]
kWindow       = 0.8          # Thermal conductivity of glass [W/mK]
U_windows = ((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm))**-1
m_windows = windowDensity * windows_area * LWindow

# Walls (concrete)
lenHouse    = 15             # House length [m]
widHouse    = 8              # House width [m]
htHouse     = 2.6            # House height [m]
LWall       = 0.25           # Wall thickness [m]
wallDensity = 2400           # Density [kg/m^3]
c_wall      = 750            # Specific heat [J/kgK]
kWall       = 0.14           # Thermal conductivity [W/mK]
walls_area = 2*(lenHouse + widHouse) * htHouse - windows_area
U_wall = ((1/h_air_wall) + (LWall /kWall) + (1/h_wall_atm))**-1
m_walls = wallDensity * walls_area * LWall

# Roof (glass fiber)
pitRoof     = 40/180/pi      # Roof pitch (40 deg)
LRoof       = 0.2            # Roof thickness [m]
roofDensity = 2440           # Density of glass fiber [kg/m^3]
c_roof      = 835            # Specific heat of glass fiber [J/kgK]
kRoof       = 0.04           # Thermal conductivity of glass fiber [W/mK]
roof_Area = 2 * (widHouse/(2*cos(pitRoof))*lenHouse)
U_roof = ((1/h_air_roof) + (LRoof/kRoof) + (1/h_roof_atm))**-1
m_roof = roofDensity * roof_Area * LRoof


m_air = airDensity * lenHouse * widHouse * htHouse

mc_T = m_air*c_air + m_roof*c_roof + m_windows*c_window + m_walls*c_wall

# PV
n_modules = 10
module_power_ref = 0.315
module_power = 0.400

################################   Simulations ################################

import matplotlib.pyplot as plt
import csvreader
from numpy import savetxt

start = time.time()    # The timer is initializad.
start_day = 0
end_day = 365

CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataP_Load = csvreader.read_data(csv='Load_Profile_15min.csv', address='', delim=',')
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataPV.data2array()
CSVDataTamb.data2array()
CSVDataP_Load.data2array()
CSVDataRad.data2array()
P_PV = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]
T_amb = [i[0]+273 for i in CSVDataTamb.ar[start_day*24*4:end_day*24*4]]
P_Load = [i for i in CSVDataP_Load.ar[0]]
a = arange(0,len(CSVDataRad.ar),15)
G = array([CSVDataRad.ar[i][0] for i in a])


# Initial conditions
t = array([0])                      # In s
dt = 60*15                          # In s
t_final = int((end_day - start_day)*24*3600/dt)         # In s

# TESS
T_TESS = array([75 + 273])          # In K
Qdot_TESS = array([0, 0])              # In W
TESS_active = False

# HP
HP_active = False
HP_status = HP_Power(HP_active)
P_HP = array([HP_status[0],HP_status[0],HP_status[0]])       # In kW
Qdot_HP = array([HP_status[1],HP_status[1],HP_status[1]])    # In kW


# Thermal demand
T_0 = 20 + 273
T_in = array([T_0, T_0, T_0])          # In K
T_set_day = [17+273]*int((6-0)*4) + [20+273]*int((22-6)*4)+ [17+273]*int((24-22)*4)
T_set = array(T_set_day*(end_day-start_day))
Qdot_Losses = array([0, 0, 0])
Qdot_SC_TESS = array([0])

# Solar Collectors
SC_active = False
Qdot_SC = array([0, 0, 0])                # In W


# BESS
SoC_BESS = array([.50])                     # In %
P_BESS = array([0])                         # In kW


# Grid
P_Grid = array([0])                         # In kW



for step in range(t_final-1):
    
    ################### Thermal carrier ######################
    
    ## Thermal losses
    # Roof
    Qdot_roof = heat_loss(U_roof, roof_Area, T_in[step], T_amb[step])
    
    # Windows
    Qdot_windows = heat_loss(U_windows, windows_area, T_in[step], T_amb[step])
    
    # Walls
    Qdot_wall = heat_loss(U_wall, walls_area, T_in[step], T_amb[step])
    
    Qdot_losses = Qdot_roof + Qdot_windows + Qdot_wall


    # Solar Collector
    if T_in[-1] < T_set[step-1]: # and Qdot_losses > Qdot_SC[-1]:  
        SC_active = True
        
    else:
        SC_active = False 

    Qdot_SC = append(Qdot_SC, Qdot_SolarCollector(SC_active, G[step]))
    
    # TESS
    if T_in[-2] < T_set[step-2]: # and Qdot_losses > Qdot_SC[-1]:         
        TESS_active = True
        
    else:
        TESS_active = False     
   
    Qdot_SC_TESS = append(Qdot_SC_TESS, Qdot_SolarCollector(not SC_active, G[step]))
    TESS_state = update_TESS(TESS_active, T_TESS[step], T_soil = T_amb[step], Qdot_SC = Qdot_SC_TESS[-1], dt = dt)
    T_TESS = append(T_TESS, TESS_state[0])
    Qdot_TESS =  append(Qdot_TESS, TESS_state[1])

    
    # HP
    if T_in[-3] < T_set[step-3]: # and Qdot_losses > Qdot_SC[-1]:
        HP_active = True
    else:
        HP_active = False
        
    HP_state = HP_Power(HP_active)
    P_HP = append(P_HP, HP_state[0])
    Qdot_HP = append(Qdot_HP, HP_state[1])        
    
    # Thermal demand
    
    # Extra code to estimate the thermal load to keep the house at the initial condition T_0
#    if new_house_Temperature(T_in[step], Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_SC = Qdot_SC[step], dt = dt) > (T_0 + 273):
#        T_in = append(T_in, new_house_Temperature(T_in[step], Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_SC = Qdot_SC[step], dt = dt))
#    else:
#        T_in = append(T_in, T_0 + 273)
    
    # Real code
    
#    print("step", step, ", Tin ", T_in[-1], Qdot_roof, Qdot_windows, Qdot_wall, Qdot_HP[-1], Qdot_TESS[-1], Qdot_SC[-1])
    T_in = append(T_in, new_house_Temperature(T_in[-1], Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_HP = Qdot_HP[-1], Qdot_TESS = Qdot_TESS[-1], Qdot_SC = Qdot_SC[-1], dt = dt))
    t = append(t, dt*step/3600)


    # Extra code to estimate the thermal load to keep the house at the initial condition T_0
#    if (Qdot_roof+Qdot_windows+Qdot_wall) > 0:
#        Qdot_Losses = append(Qdot_Losses, 0)
#    else:
#        Qdot_Losses = append(Qdot_Losses, -(Qdot_roof+Qdot_windows+Qdot_wall)) 
    
    # Real code    
    Qdot_Losses = append(Qdot_Losses, -(Qdot_roof+Qdot_windows+Qdot_wall))


    ################### Electric carrier ######################

    

    # BESS
    [SoC_BESS_state, P_Grid_state, P_BESS_state] = update_BESS(SoC_BESS[-1], P_Load[step] + P_HP[-1], P_PV[step])
    P_BESS = append(P_BESS, P_BESS_state)    
    SoC_BESS = append(SoC_BESS, SoC_BESS_state)
    P_Grid = append(P_Grid, P_Grid_state)


end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.
###################################   CSVs  ###################################

#savetxt('T_in_15min.csv', T_in, delimiter =", ", fmt ='% s')
#savetxt('Thermal_load_15min.csv', Qdot_Losses, delimiter =", ", fmt ='% s')


##################################   Plots  ###################################

plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman"
})


plt.figure(1)
plt.plot(t, [i-273 for i in T_TESS])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature in the TESS, $T_{TESS}$, [°C]')
#plt.title('Temperature in the TESS')
plt.show()

plt.figure(2)
plt.plot(t, [i/1000 for i in Qdot_TESS[1:]])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the TESS, $\dot{Q}_{TESS}$, [kW]')
#plt.title('Thermal power provided by the TESS')
plt.show()

plt.figure(3)
plt.plot(t, [i-273 for i in T_amb], 'r', label='Ambient temperature, $T_{amb}$')
plt.plot(t, [i-273 for i in T_in[2:]], 'b', label='Temperature inside the house, $T_{in}$')
plt.plot(t, [i-273 for i in T_set], 'g', label='Setpoint temperature inside the house, $T_{set}$')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature [°C]')
#plt.title('Temperature inside the house')
plt.show()    
#
#plt.figure(4)
#plt.plot(t, [i/1000 for i in Qdot_Losses[2:]])
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylim([0, 1.8])
#plt.ylabel('Thermal losses, $\dot{Q}_{L}$, [kW]')
##plt.title('Thermal losses in the house')
#plt.show()    
#
plt.figure(5)
plt.plot(t, [i/1000 for i in Qdot_HP[2:]])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the HP, $\dot{Q}_{HP}$, [kW]')
#plt.title('Heat provided by the HP')
plt.show()

plt.figure(6)
plt.plot(t, [i/1000 for i in Qdot_SC[2:]], 'b', label='To the thermal load')
plt.plot(t, [i/1000 for i in Qdot_SC_TESS], 'r', label='To the TESS')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the SC, $\dot{Q}_{SC}$, [kW]')
#plt.title('Thermal power provided by the SC')
plt.show()
#
plt.figure(7)
plt.plot(t, P_Load)
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Electric load, $P_{L}$, [kW]')
#plt.title('Electric load')
plt.show()  
#
#plt.figure(8)
#plt.plot(t, [i/1000 for i in P_HP[2:]])
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylabel('Power of the HP, $P_{HP}$, [kW]')
##plt.title('Power consumed by the HP')
#plt.show()

plt.figure(9)
plt.plot(t, P_PV[:-1])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('PV power load, $P_{PV}$, [kW]')
#plt.title('Electric load')
plt.show()  

plt.figure(10)
plt.plot(t, SoC_BESS)
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('State-of-charge of the BESS, $SoC_{BESS}$ [%]')
#plt.title('SoC of the BESS')
plt.show()

plt.figure(11)
plt.plot(t, P_BESS, 'b', label='Power delivered by the BESS')
#plt.legend(loc='center right')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Power of the BESS, $P_{BESS}$, [kW]')
#plt.title('Power of the BESS')
plt.show()

plt.figure(12)
plt.plot(t, P_Grid, 'b', label='Power consumed from the grid')
#plt.legend(loc='center right')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Power from the grid, $P_{G}$, [kW]')
#plt.title('Power from the grid')
plt.show()