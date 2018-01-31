# Size of variable arrays:
sizeAlgebraic = 12
sizeStates = 4
sizeConstants = 5
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "Cm in component membrane (microF)"
    legend_algebraic[4] = "i_Na in component sodium_channel (nanoA)"
    legend_algebraic[10] = "i_K in component potassium_channel (nanoA)"
    legend_algebraic[11] = "i_Leak in component leakage_current (nanoA)"
    legend_constants[1] = "g_Na_max in component sodium_channel (microS)"
    legend_algebraic[0] = "g_Na in component sodium_channel (microS)"
    legend_constants[2] = "E_Na in component sodium_channel (millivolt)"
    legend_states[1] = "m in component sodium_channel_m_gate (dimensionless)"
    legend_states[2] = "h in component sodium_channel_h_gate (dimensionless)"
    legend_algebraic[1] = "alpha_m in component sodium_channel_m_gate (per_second)"
    legend_algebraic[5] = "beta_m in component sodium_channel_m_gate (per_second)"
    legend_algebraic[2] = "alpha_h in component sodium_channel_h_gate (per_second)"
    legend_algebraic[6] = "beta_h in component sodium_channel_h_gate (per_second)"
    legend_algebraic[8] = "g_K1 in component potassium_channel (microS)"
    legend_algebraic[9] = "g_K2 in component potassium_channel (microS)"
    legend_states[3] = "n in component potassium_channel_n_gate (dimensionless)"
    legend_algebraic[3] = "alpha_n in component potassium_channel_n_gate (per_second)"
    legend_algebraic[7] = "beta_n in component potassium_channel_n_gate (per_second)"
    legend_constants[3] = "g_L in component leakage_current (microS)"
    legend_constants[4] = "E_L in component leakage_current (millivolt)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt m in component sodium_channel_m_gate (dimensionless)"
    legend_rates[2] = "d/dt h in component sodium_channel_h_gate (dimensionless)"
    legend_rates[3] = "d/dt n in component potassium_channel_n_gate (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -87
    constants[0] = 12
    constants[1] = 400000
    constants[2] = 40
    states[1] = 0.01
    states[2] = 0.8
    states[3] = 0.01
    constants[3] = 75
    constants[4] = -60
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = (100.000*(-states[0]-48.0000))/(exp((-states[0]-48.0000)/15.0000)-1.00000)
    algebraic[5] = (120.000*(states[0]+8.00000))/(exp((states[0]+8.00000)/5.00000)-1.00000)
    rates[1] = algebraic[1]*(1.00000-states[1])-algebraic[5]*states[1]
    algebraic[2] = 170.000*exp((-states[0]-90.0000)/20.0000)
    algebraic[6] = 1000.00/(1.00000+exp((-states[0]-42.0000)/10.0000))
    rates[2] = algebraic[2]*(1.00000-states[2])-algebraic[6]*states[2]
    algebraic[3] = (0.100000*(-states[0]-50.0000))/(exp((-states[0]-50.0000)/10.0000)-1.00000)
    algebraic[7] = 2.00000*exp((-states[0]-90.0000)/80.0000)
    rates[3] = algebraic[3]*(1.00000-states[3])-algebraic[7]*states[3]
    algebraic[0] = (power(states[1], 3.00000))*states[2]*constants[1]
    algebraic[4] = (algebraic[0]+140.000)*(states[0]-constants[2])
    algebraic[8] = 1200.00*exp((-states[0]-90.0000)/50.0000)+15.0000*exp((states[0]+90.0000)/60.0000)
    algebraic[9] = 1200.00*(power(states[3], 4.00000))
    algebraic[10] = (algebraic[8]+algebraic[9])*(states[0]+100.000)
    algebraic[11] = constants[3]*(states[0]-constants[4])
    rates[0] = -(algebraic[4]+algebraic[10]+algebraic[11])/constants[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = (100.000*(-states[0]-48.0000))/(exp((-states[0]-48.0000)/15.0000)-1.00000)
    algebraic[5] = (120.000*(states[0]+8.00000))/(exp((states[0]+8.00000)/5.00000)-1.00000)
    algebraic[2] = 170.000*exp((-states[0]-90.0000)/20.0000)
    algebraic[6] = 1000.00/(1.00000+exp((-states[0]-42.0000)/10.0000))
    algebraic[3] = (0.100000*(-states[0]-50.0000))/(exp((-states[0]-50.0000)/10.0000)-1.00000)
    algebraic[7] = 2.00000*exp((-states[0]-90.0000)/80.0000)
    algebraic[0] = (power(states[1], 3.00000))*states[2]*constants[1]
    algebraic[4] = (algebraic[0]+140.000)*(states[0]-constants[2])
    algebraic[8] = 1200.00*exp((-states[0]-90.0000)/50.0000)+15.0000*exp((states[0]+90.0000)/60.0000)
    algebraic[9] = 1200.00*(power(states[3], 4.00000))
    algebraic[10] = (algebraic[8]+algebraic[9])*(states[0]+100.000)
    algebraic[11] = constants[3]*(states[0]-constants[4])
    return algebraic

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0,2,1000)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    #pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.plot(voi,states[0])
    pylab.xlabel(legend_voi)
    #pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)