# BNCE Docs

**Written by:** Vinay Williams, vinaywilliams00@gmail.com

This program uses the Rao Method to compute the nozzle's contour, there is plans to add functionality for its determination using the method of characteristics

References:

* Rao, G.V.R., 1958. Exhaust nozzle contour for optimum thrust. *Journal of Jet Propulsion*, *28*(6), pp.377-382.
* Davis, K., Fortner, E., Heard, M., McCallum, H. and Putzke, H., 2015. Experimental and computational investigation of a dual-bell nozzle. In *53rd AIAA Aerospace Sciences Meeting* (p. 0377).

### Current Functionality

* Determines the contour of a bell nozzle
* Determines the converging region angle and contour nearest the throat
* Plots the contour

### Future Functionality 

* Addition of Performance Statistics
* Additions of a secondary solver, i.e Method of Characteristics
* Options for use of .xls file for supply of variables

### Supplied Variables

```
mdot: Mass Flow Rate, Default Value = 0
G: Specific Heat Ratio
theta_1st: Initial Angle for Coverging Curve from throat
theta_nth: Initial Angle for Diverging Curve form throat
RT: Radius of Throat
TC: Chamber Temperature
PA: Ambient Pressure
PC: Chamber Pressure
k: Useful Length Fraction
```

### Important Calculated Values

``` 
TRatio: Temperature Ratio
G1: Constant = G+1
G2: Constant = G-1
PR_Throat: Pressure Ratio
TR_throat: Temperature Ratio at throat
T: Throat Temperature
P: Throat Pressure
AT: Throat Area
e: Expansion ratio
M: Exit Mach Number
Ae: Exit Area
Re: Exit Radius
L: Length of Nozzle

y: y points for nozzle contour
x: x points for nozzle contour
```

