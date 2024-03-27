Ray-tracing software (written in matlab) used to simulate the effect of reabsorption of emission light by trapped charges in cylindrical SrAl2O4:Eu,Dy single crystals.

**The following sample-specific properties can be changed:**
- the degree of reabsorption (effectStrength),
- the size (length L, radius R) of the crystal
- refractive index of the crystal (n1)
- the refractive index of the surrounding medium (n2)

**The following properties of the trapping process can be changed:** 
- trap depth (ET)
- frequency factor (sET)
- trap density (NT)

**The flow of the simulation is as follows:**
- Determine how many charges are detrapped based on trap concentration, trap depth and frequency factor.
- For every detrapped charge:
        * create new photon at random location, in random direction.
        * create random absorption path length from distribution based on concentration and effectStrength
        * calculate intersection between direction of photon and crystal surface.
        * calculate path length
        * if path length exceeded, determine position of reabsorption, decrease trapped charges by 1 and restart
        * if not calculate angle of incidence and check if internal reflection or not
        * if internal reflection => calculate new direction photon, if not => photon escapes, intensty at this timestep augmented by 1
- Store amount of escaped photons, move to next timestep until number of trapped charges is below threshold value. 

     
