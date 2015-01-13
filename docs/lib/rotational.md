



# Rotational mechanics

Library to model 1-dimensional, rotational mechanical systems

Rotational provides 1-dimensional, rotational mechanical components to
model in a convenient way drive trains with frictional losses.

These components are modeled after the Modelica.Mechanics.Rotational
library.

NOTE: these need more testing.




## Inertia

1D-rotational component with inertia

Rotational component with inertia and two rigidly connected flanges. 

```julia
Inertia(flange_a::Flange, flange_b::Flange, J::Real)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `J::Real` : Moment of inertia [kg.m^2]


[Sims/src/../lib/rotational.jl:42](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L42)



## Disc

1-dim. rotational rigid component without inertia, where right flange is rotated by a fixed angle with respect to left flange

Rotational component with two rigidly connected flanges without
inertia. The right flange is rotated by the fixed angle "deltaPhi"
with respect to the left flange.

```julia
Disc(flange_a::Flange, flange_b::Flange, deltaPhi)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `deltaPhi::Signal` : rotation of left flange with respect to right flange (= flange_b - flange_a) [rad]


[Sims/src/../lib/rotational.jl:77](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L77)



## Spring

Linear 1D rotational spring

A linear 1D rotational spring. The component can be connected either
between two inertias/gears to describe the shaft elasticity, or
between a inertia/gear and the housing (component Fixed), to describe
a coupling of the element with the housing via a spring.

```julia
Spring(flange_a::Flange, flange_b::Flange, c::Real, phi_rel0 = 0.0)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `c`: spring constant [N.m/rad]
* `phi_rel0` : unstretched spring angle [rad]


[Sims/src/../lib/rotational.jl:105](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L105)



## Damper

Linear 1D rotational damper

Linear, velocity dependent damper element. It can be either connected
between an inertia or gear and the housing (component Fixed), or
between two inertia/gear elements.

```julia
Damper(flange_a::Flange, flange_b::Flange, d::Signal)
Damper(flange_a::Flange, flange_b::Flange, hp::HeatPort, d::Signal)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `hp::HeatPort` : heat port [K]
* `d`: 	damping constant [N.m.s/rad]


[Sims/src/../lib/rotational.jl:136](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L136)



## SpringDamper

Linear 1D rotational spring and damper in parallel

A spring and damper element connected in parallel. The component can
be connected either between two inertias/gears to describe the shaft
elasticity and damping, or between an inertia/gear and the housing
(component Fixed), to describe a coupling of the element with the
housing via a spring/damper.

```julia
SpringDamper(flange_a::Flange, flange_b::Flange, c::Signal, d::Signal)
SpringDamper(flange_a::Flange, flange_b::Flange, hp::HeatPort, c::Signal, d::Signal)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `hp::HeatPort` : heat port [K]
* `c`: 	spring constant [N.m/rad]
* `d`: 	damping constant [N.m.s/rad]


[Sims/src/../lib/rotational.jl:172](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L172)



## IdealGear

Ideal gear without inertia

This element characterices any type of gear box which is fixed in the
ground and which has one driving shaft and one driven shaft. The gear
is ideal, i.e., it does not have inertia, elasticity, damping or
backlash. If these effects have to be considered, the gear has to be
connected to other elements in an appropriate way.

```julia
IdealGear(flange_a::Flange, flange_b::Flange, ratio)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `ratio` : transmission ratio (flange_a / flange_b)


[Sims/src/../lib/rotational.jl:293](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L293)




# Miscellaneous




## MBranchHeatPort

Wrap argument `model` with a heat port that captures the power
generated by the device. This is vectorizable.

```julia
MBranchHeatPort(flange_a::Flange, flange_b::Flange, hp::HeatPort,
                model::Function, args...)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `hp::HeatPort` : Heat port [K]                
* `model::Function` : Model to wrap
* `args...` : Arguments passed to `model`  


[Sims/src/../lib/rotational.jl:331](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L331)




# Sensors




## SpeedSensor

Ideal sensor to measure the absolute flange angular velocity

Measures the absolute angular velocity w of a flange in an ideal way
and provides the result as output signal w.

```julia
SpeedSensor(flange::Flange, w::Signal)
```

### Arguments

* `flange::Flange` : left flange of shaft [rad]
* `w::Signal`: 	absolute angular velocity of the flange [rad/sec]


[Sims/src/../lib/rotational.jl:377](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L377)



## AccSensor

Ideal sensor to measure the absolute flange angular acceleration

Measures the absolute angular velocity a of a flange in an ideal way
and provides the result as output signal a.

```julia
SpeedSensor(flange::Flange, a::Signal)
```

### Arguments

* `flange::Flange` : left flange of shaft [rad]
* `a::Signal`: 	absolute angular acceleration of the flange [rad/sec^2]


[Sims/src/../lib/rotational.jl:400](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L400)




# Sources




## SignalTorque

Input signal acting as external torque on a flange

The input signal tau defines an external torque in [Nm] which acts
(with negative sign) at a flange connector, i.e., the component
connected to this flange is driven by torque tau.

```julia
SignalTorque(flange_a::Flange, flange_b::Flange, tau::Signal)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `tau` : Accelerating torque acting at flange_a relative to flange_b
  (normally a support); a positive value accelerates flange_a


[Sims/src/../lib/rotational.jl:437](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L437)



## QuadraticSpeedDependentTorque

Quadratic dependency of torque versus speed

Model of torque, quadratic dependent on angular velocity of flange.
Parameter TorqueDirection chooses whether direction of torque is the
same in both directions of rotation or not.

```julia
QuadraticSpeedDependentTorque(flange_a::Flange, flange_b::Flange,
                              tau_nominal::Signal, TorqueDirection::Bool, w_nominal::Signal)
```

### Arguments

* `flange_a::Flange` : left flange of shaft [rad]
* `flange_b::Flange` : right flange of shaft [rad]
* `tau_nominal::Signal` : nominal torque (if negative, torque is acting as a load) [N.m]
* `TorqueDirection::Bool` : same direction of torque in both directions of rotation
* `AngularVelocity::Signal` : nominal speed [rad/sec]


[Sims/src/../lib/rotational.jl:465](https://github.com/tshort/Sims.jl/tree/d39a15c1969c6fad87a4a7ab7f25088963690512/src/../lib/rotational.jl#L465)
