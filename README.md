# cryoH2vehicle simulation 

Using cryogenic H2 onboard a vehicle has the proven advantage to dramatically increase the density of the fluid as compared to room temperature. Under certain par/drive/fill scenarios, however, the pressure in vessel can reach its rated value and the hydrogen may start leaking through the pressure relief device. 

A FORTRAN code was written to simulate the variations of the thermodynamics state of the hydrogen inside the vessel (pressure, temperature, density, para/ortho) as a function of the vessel design (volume, rated pressure, aspect ratio, heat transfer – including the influence of outside temperature) and the duty cycle the vehicle experiences, in terms of parking, driving, and refilling (assuming the vehicle always refills when it reaches a fixed amount, using a dispensing system of fixed outlet entropy or temperature). As a result, large parametric studies can be performed (influence of outside temperature, volume, rated pressure, fill pressure, dispensing entropy…) for various duty cycles, to be specified by the user; and typical states of charge and boil-off losses can be estimated. The code relies on the conservation of mass and energy, and includes real gas equations of state, para/ortho kinetics, temperature dependent heat capacity for the materials, and a model for tank design.

# How to use

See [USER MANUAL](https://github.com/LLNL/cryoH2vehicle/blob/master/UserManual.pdf)

# License

The software is released under the LLNL LGPG license. For more details see the [LICENSE](https://github.com/LLNL/cryoH2vehicle/blob/master/LICENSE.md/LICENSE.md) file.

LLNL-CODE-750958
