import numpy as np
from Typing import Tuple

def power(velocity, 
          mass, 
          specific_energy=0.98, 
          drag_coefficient=0.24, 
          air_density=1.225, 
          windspeed=0, 
          hill_gradient=0) -> Tuple[float, float, float, float]:
    """
    Running power needed to achieve given velocity with given mass, etc.
    
    :param float velocity: The running speed [m/s]
    :param float mass: The total weight of the runner [kg]
    :param float specific_energy: The specific cost of running, default 0.98 [kJ/kg/km]
    :param float drag_coefficient: The air drag coefficient cdA, default 0.25 [m²]
    :param float air_density: The air density, default standard conditions 1.225 [kg/m³]
    :param float windspeed: The head wind speed, default 0 [m/s]
    :param float hill_gradient: The gradient of the hill, default 0 [%]
    :return: Total power, running resistance power, air resistance power, climbing power
    :rtype: Tuple[float, float, float, float]
    """        
    g = 9.81 # [m/s²] 
    eta = (45.6 + 1.1622*hill_gradient)/100 # Minetti formula, book p. 88
    
    power_running  = specific_energy*mass*velocity
    power_air      = 0.5*air_density*drag_coefficient*(velocity + windspeed)**2*velocity 
    power_climbing = (hill_gradient/100)*mass*g*velocity*eta
    
    power = power_running + power_air + power_climbing
    
    return power, power_running, power_air, power_climbing

def velocity(power, 
          mass, 
          specific_energy=0.98, 
          drag_coefficient=0.24, 
          air_density=1.225, 
          windspeed=0, 
          hill_gradient=0) -> float:
    """
    Running speed given power, mass, etc.
    
    :param float velocity: The running speed [m/s]
    :param float mass: The total weight of the runner [kg]
    :param float specific_energy: The specific cost of running, default 0.98 [kJ/kg/km]
    :param float drag_coefficient: The air drag coefficient cdA, default 0.25 [m²]
    :param float air_density: The air density, default standard conditions 1.225 [kg/m³]
    :param float windspeed: The head wind speed, default 0 [m/s]
    :param float hill_gradient: The gradient of the hill, default 0 [%]
    :return: Running velocity
    :rtype: float
    """    
    
    g = 9.81 # [m/s²] 
    eta = (45.6 + 1.1622*hill_gradient)/100 # Minetti formula, book p. 88
    

    # https://en.wikipedia.org/wiki/Cubic_equation#General_cubic_formula
    a = 0.5*air_density*drag_coefficient
    b = a*2*windspeed
    c = a*windspeed**2 + eta*hill_gradient*mass*g + specific_energy*mass
    d = -power
    
    delta_0 = b**2 - 3*a*c
    delta_1 = 2*b**3 - 9*a*b*c + 27*a**2*d
    
    C = np.cbrt((delta_1 + np.sqrt(delta_1**2 - 4*delta_0**3))/2)
    
    velocity = -1/(3*a)*(b + C + delta_0/C)
    
    return velocity