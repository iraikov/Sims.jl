# Industrial Control Systems package
#
#
# Copyright 2012 Electric Power Research Institute
# 
# MIT License - this license change allowed by the Modelica License
#
# Derived from the following:
#
# Modelica library for Industrial Control Systems
#  
# Licensed by Marco Bonvini and Alberto Leva under the Modelica License 2 Copyright © 2010-2012
# 
# Version: 1.0.0
# Date: May, 15, 2012
# 
# Marco Bonvini <bonvini@elet.polimi.it, bonvini.m@gmail.com>
# Alberto Leva <leva@elet.polimi.it>
#
# More info and download:
#    http://marcobonvini.altervista.org/ics.php
# Online docs:
#    http://marcobonvini.altervista.org/IndustrialControlSystems/help/IndustrialControlSystems.html
# Paper:
#    Marco Bonvini and Alberto Leva, "A Modelica Library for
#    Industrial Control Systems", Proceedings of the 9th International
#    Modelica Conference, September, 2012
#    http://www.ep.liu.se/ecp/076/048/ecp12076048.pdf 
#

require("Sims")
require("Options")

module ICS       # IndustrialControlSystems


using Sims
using OptionsMod

export PID, FirstOrder, ComplexPoles, Step

# Discrete integrator
function fIntegrator(alfa::Float64,  # Parametro di discretizzazione
                     u::Float64,     # Current input
                     u_pre::Float64, # Previous input
                     y_pre::Float64, # Previous output
                     Ts::Float64,    # Sampling time [s]
                     k::Float64)     # Gain
    A = alfa*Ts
    B = Ts - alfa*Ts
    y = A*k*u + B*k*u_pre + y_pre
    return y
end

# Discrete first order process : k/(1+s*tau)
function f1Pole(alfa::Float64,  # Parametro di discretizzazione
                u::Float64,     # Current input
                u_pre::Float64, # Previous input
                y_pre::Float64, # Previous output
                Ts::Float64,    # Sampling time [s]
                k::Float64,     # Gain
                tau::Float64)   # Pole
    A = alfa*Ts
    B = Ts - alfa*Ts
    C = A + tau
    D = B - tau
    y = (A*k*u + B*k*u_pre - D*y_pre)/C
    return y
end 

pre(d::Discrete) = d.pre 

function PID(SP::Signal, PV::Signal, CS::Signal, opts::Options)
    @defaults opts begin 
        # common entries
        Ts = 0.0  # Sampling time (if <= 0 continuous time)
        method = "BE" # Discretisation method
        AntiWindup = false # Flag that enables the antiwindup feature
        CSmin = 0.0    # Minimum value of the CS
        CSmax = 1.0    # maximum value of the CS
        CS_start = 0.0 # output initial value
        eps = 1e-6 # small time constant that represents the time for switching between auto and tracking mode
        # specific to PID
        tr = nothing   # track reference signal
        bias = 0.0     # biasing signal
        Kp = 5.0  # Proportional gain
        Ti = 1.0  # Integral time
        Td = 1.0  # Derivative time
        Tt = max(1e3*eps,sqrt(Ti*Td))  # Reset time for the integrator (it should be larger than Td and smaller than Ti)
        N = 10.0  # Derivative filter ratio
        b =  1.0  # Set point weighting for the proportional action
        c =  1.0  # Set point weighting for the derivative action
    end
    ts = tr != nothing  # track switch signal

    alfa = ["BE" => 1.0, "FE" => 0.0, "TU" => 0.5][method]
    
    assert(Ti>0, "Integral time `Ti` must be positive")
    assert(Tt>0, "Reset time `Tt` must be positive")
    assert(Td>=0,"Derivative time `Td` must be >= 0")
    
    if Ts == 0.0 # Continuous version
        Paction = Unknown() # proportional action
        Iaction = Unknown() # integral action
        Daction = Unknown() # derivative action
        Fder    = Unknown() # filtered weighted error
        satin   = Unknown() # saturation input
        satDiff = Unknown() # difference between saturation input and output signals
        cs      = Unknown() # Control Signal before the biasing

        # EQUATIONS
        {
         Paction - Kp*(b*SP - PV)

         # integral action
         if !ts
           -der(Iaction) + Kp/Ti*(SP - PV) + 1/Tt*satDiff
         else
           Iaction + eps*der(Iaction) - cs
         end

         # derivative action
         -Daction + Kp*N*(c*SP - PV) - Fder
         -Fder    + ((Td > 0.0) ? (-Td/N*der(Fder) + Kp*N*(c*SP - PV)) : 0.0)

         # sum of the three actions
         -satin + Iaction + Paction + (Td>0 ? Daction : 0.0)

         # Anti-windup
         -cs + (AntiWindup ? max(CSmin,min(CSmax,ts ? tr : satin)) : (ts ? tr : satin))

         # Signal for preventing antiwindup (integrator reset signal)
         -satDiff + cs - satin

         # biasing the output
         -CS + (AntiWindup ? max(CSmin,min(CSmax,cs + bias)) : (cs + bias))
         ## NOTE: May need a Limiter instead of the max and min functions above.
         ##       min/max don't trigger events.
        }
    else   # Discrete version
        Paction = Discrete() # proportional action
        Iaction = Discrete() # integral action
        Daction = Discrete() # derivative action
        Fder    = Discrete() # filtered weighted error
        satin   = Discrete() # saturation input
        satDiff = Discrete() # difference between saturation input and output signals
        cs      = Discrete() # Control Signal before the biasing
        SPd     = Discrete() 
        PVd     = Discrete() 
        {
         Event(sin(MTime / Ts * 2pi), {
             reinit(Paction, Kp*(b*SP - PV))
             # Signal for preventing antiwindup (integrator reset signal)
             reinit(satDiff, cs - satin)
             reinit(satin, Iaction + Paction + (Td>0 ? Daction : 0))
             # Anti-windup
             reinit(cs, AntiWindup ? max(CSmin,min(CSmax,(ts ? tr : satin))) : (ts ? tr : satin))
             reinit(Iaction, !ts ?
                    fIntegrator(alfa, Kp/Ti*(SP - PV) + 1/Tt*satDiff, Kp/Ti*(SPd.pre - PVd.pre) + 1/Tt*satDiff.pre, Iaction, Ts, 1) :
                    f1Pole(alfa, cs, cs.pre, Iaction.pre, Ts, 1, eps))
             reinit(Fder, Td > 0 ?
                    f1Pole(alfa, c*SP - PV, c*SPd.pre - PVd.pre, Fder.pre, Ts, Kp*N, Td/N) :
                    0)
             reinit(Daction, Kp*N*(c*SP - PV) - Fder)
             reinit(SPd, SP)
             reinit(PVd, PV)
         })
         # biasing the output
         CS - (AntiWindup ? max(CSmin,min(CSmax,cs + pre(bias))) : (cs + pre(bias))) 
        }
    end
end

function FirstOrder(u::Signal, y::Signal, opts::Options)
    @defaults opts begin 
        tau = 2.0  # pole's time constant
        mu = 1.0   # Gain
        y_start = 0.0 # output initial value
    end
    y.value = y_start
    {
     y + tau*der(y) - mu*u
    }
end

function ComplexPoles(u::Signal, y::Signal, opts::Options)
    @defaults opts begin 
        xi = 0.8       #  Damping coefficient
        omegan = 0.1   #  Natural frequency
        mu = 1.0       #  Gain
        y_start = 0.0  #  Output initial value
        dy_start = 0.0 #  Slope initial value
    end
    assert(xi >= 0 && xi <= 1,"Dumping coefficient `xi` must be between 0 and 1.")
    dy = Unknown(dy_start)
    y.value = y_start
    {
     dy - der(y)
     -mu*u + y + 2*xi/omegan*dy + der(dy)/(omegan^2)
    }
end

## function TransferFunction(u::Signal, y::Signal, opts::Options)
##     @defaults opts begin 
##         num = [1.0,1]   # Numerators coefficients (4*s + 2) is [4,2]
##         den = [2,1,1]   # Denumerators coefficients (4*s + 2) is [4,2]
##         y_start = 0.0  #  Output initial value
##     end
##     assert(xi >= 0 && xi <= 1,"Dumping coefficient `xi` must be between 0 and 1.")
##     dy = Unknown(dy_start)
##     y.value = y_start
##     {
##      dy - der(y)
##      -mu*u + y + 2*xi/omegan*dy + der(dy)/(omegan^2)
##     }
## end

function Step(y::Signal, opts::Options)
    @defaults opts begin
        height = 1.0  
        offset = 0.0 
        startTime = 0.0
    end
    ymag = Discrete(offset)
    {
     y - ymag  
     Event(MTime - startTime,
           {reinit(ymag, offset + height)},   # positive crossing
           {reinit(ymag, offset)})            # negative crossing
    }
end
Step(y::Signal) = Step(y, Options())

end # module ICS


using ICS
using Sims
using OptionsMod

function ex_ProcessControl()
    step_u = Unknown("step_u")
    FO_u = Unknown("FO_u")
    FO_y = Unknown("FO_y")
    FO_uncontrolled = Unknown("FO_uncontrolled")
    Osc_u = Unknown("Osc_u")
    Osc_y = Unknown("Osc_y")
    Osc_uncontrolled = Unknown("Osc_uncontrolled")
    {
     Step(step_u, @options(startTime => 1.0))
     FirstOrder(step_u, FO_uncontrolled, @options(tau => 5.0, mu => 1.0))
     FirstOrder(FO_u, FO_y, @options(tau => 5.0, mu => 1.0))
     PID(step_u, FO_y, FO_u,
         @options(AntiWindup => true,
         CSmin => 0.0,
         CSmax => 5.0,
         Kp => 10.0,
         Td => 0.8,
         N => 8.0,
         Ti => 3.0))
     ComplexPoles(step_u, Osc_uncontrolled, @options(xi => 0.3, omegan => 0.5, mu => 1.0))
     ComplexPoles(Osc_u, Osc_y, @options(xi => 0.3, omegan => 0.5, mu => 1.0))
     PID(step_u, Osc_y, Osc_u,
         @options(AntiWindup => true,
         CSmin => 0.0,
         CSmax => 5.0,
         Kp => 10.0,
         Td => 1.0,
         Ti => 2.5))
    }
end

function sim_ProcessControl()
    m = ex_ProcessControl()
    f = elaborate(m)
    s = create_sim(f)
    y = sim(s, 40.0)
    wplot(y, "ProcessControl.pdf")
end

function ex_PI_ContinuousVsDigital()
    step_u = Unknown("step_u")
    disturb = Unknown("disturb")
    cont_u = Unknown("cont_u")
    cont_y = Unknown("cont_y")
    cont_cs = Unknown("cont_cs")
    disc_u = Unknown("disc_u")
    disc_y = Unknown("disc_y")
    disc_cs = Unknown("disc_cs")
    {
     Step(step_u,  @options(startTime => 5.0, height => 0.5))
     Step(disturb, @options(startTime => 30.0, height => -0.2))
     ## TransferFunction(cont_u, cont_y, @options(num => [12.0,1], den => [20.0, 12, 1]))
     ## TransferFunction(disc_u, disc_y, @options(num => [12.0,1], den => [20.0, 12, 1]))
     TransferFunction(cont_u, cont_y, [12.0,1], [20.0, 12, 1])
     TransferFunction(disc_u, disc_y, [12.0,1], [20.0, 12, 1])
     disturb + cont_cs - cont_u
     disturb + disc_cs - disc_u
     PID(step_u, cont_y, cont_cs,
         @options(AntiWindup => true,
         CS_start => 0.0,
         CSmin => 0.0,
         CSmax => 2.0,
         b => 1.0,
         c => 1.0,
         Kp => 10.0,
         Td => 0.5,
         N => 8.0,
         Ti => 5.0))
     PID(step_u, disc_y, disc_cs,
         @options(AntiWindup => true,
         Ts => 0.01,
         CS_start => 0.0,
         CSmin => 0.0,
         CSmax => 2.0,
         b => 1.0,
         c => 1.0,
         Kp => 10.0,
         Td => 0.5,
         N => 8.0,
         Ti => 5.0))
    }
end

function sim_PI_ContinuousVsDigital()

    m = ex_PI_ContinuousVsDigital()
    f = elaborate(m)
    s = create_sim(f)
    y = sim(s, 40.0)
    wplot(y, "PI_ContinuousVsDigital.pdf")
    
end

function ex_TransferFunction()
    u = Unknown("u")
    y = Unknown("y")
    {
     Step(u,  @options(startTime => 5.0, height => 0.5))
     TransferFunction(u, y, [12.0,1], [20.0, 12, 1])
    }
end

function sim_TransferFunction()
    m = ex_TransferFunction()
    f = elaborate(m)
    s = create_sim(f)
    y = sim(s, 40.0)
    wplot(y, "TransferFunction.pdf")
end
