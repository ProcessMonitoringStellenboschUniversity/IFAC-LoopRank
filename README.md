# IFAC-LoopRank
Repository for code relating to submission to IFAC MMM 2019:
McCoy, JT & Auret, L. 2019. _Controller maintenance prioritisation with LoopRank for a milling and flotation plant_. Submitted.

## Overview
This paper and code provide maintenance prioritisation for control loop systems, using the LoopRank algorithm, as described in (Farenzena and Trierweiler, 2009).

The code consists of two parts: a tank network simulation, which generates synthetic data, and code related to the LoopRank algorithm, which uses MATLAB's ```partialcorr```, and Granger causality, as described in (Lindner et al., 2017a and 2017b, Wakefield et al., 2018).

Running any of ```mainFourTanks.m```, ```mainNTanks.m``` or ```mainNTanksReplicates.m``` will generate data for a simulated tank network, and evaluate the loop importances using LoopRank.

## Requirements
The code was written in MATLAB R2018a, and should not require any special toolboxes to run.

## LoopRank
LoopRank is based on PageRank, Google's algorithm for ranking web pages. It takes as input a connectivity matrix, and outputs importance scores for each page/control loop.

Here, connectivity matrices are generated using partial correlation and Granger causality, to evaluate the impact of including lagged data in causal methods, compared to the relatively simple partial correlation.

## Tank network simulation
In order to provide a ground truth for tests of the LoopRank algorithm, synthetic data were generated for a simulated tank network. Each tank has level control (with outlet flow as a manipulated variable), and some tanks also have underflow streams to other tanks in the network, with flow rate determined by tank level and an underflow constant which is proportional to tank volume. Normally-distributed noise (with a standard deviation of 0.05) is added to all measurements to represent sensor inaccuracy. The feed to the tank network fluctuates around the steady state value, representing process disturbances, calculated as an autoregressive function, as shown in equation (1). The feed flowrate at time ```k```, ```F_(feed,k)```, is a function of the flowrate at the previous time step, ```F_(feed,k-1)```, and the steady state value for the feed flowrate, ```F_(feed,SS)```, and normally distributed noise is added at each time step.

```F_(feed,k)=0.97F_(feed,k-1)+0.03F_(feed,SS)+N(0,0.15)```	(1)

The user can specify the number of tanks and the connections between them, as well as tank volumes and level setpoints, and the steady state feed flow rate to the network; example code is provided for a simple network of four tanks (mainFourTanks.m).
Example code is also provided that automatically generates a network of N tanks (the user specifies ```N```) in two banks, with level controlled flows between the tanks in the first bank, underflows from the first bank to the second bank, random level controlled flows between some of the tanks in the second bank, and random underflow recycles from the second bank to the first tank. For more detail, see the ```mainNTanks.m``` file.

### Tank model
The tank model is based on a simple volume balance (assuming incompressible fluid) over the in- and outflows from a vessel, solved using the finite difference method, as shown in equation (2). The tank level at timestep ```k+1```, ```L_(k+1)```, is calculated in terms of the tank volume, ```V```, and the values at timestep ```k``` for the feed flow rate, ```F_(in,k)```, product flow rate, ```F_(product,k)```, and underflow flow rate, ```F_(underflow,k)```.

```L_(k+1)=L_k+(F_(in,k)-F_(product,k)-F_(underflow,k))/V```	(2)

The tank model function was implemented in MATLAB, and includes logic to ensure that the tank level does not exceed 100% (product flow rate is increased to ensure this maximum level) or drop below a minimum level of 5% (product flow is set to zero and underflow calculated to maintain minimum level). For further detail, see the ```TankModel.m``` file.

### Level control
Level controllers were implemented as simple proportional integral controllers, with the controller tuning (proportional and integral gains, ```K_c``` and ```τ_i```) determined using direct synthesis (Chen and Seborg, 2002), summarised in equations (3). For each tank, the process residence time, τ,  is based on the feed flow rate, ```F_in```, and the tank volume, ```V```. The process gain, ```K```, is based on the product flow rate, ```F_product```, and the tank volume, ```V```. The tunable parameter is ```τ_c```; smaller values of the parameter result in more aggressive controllers.

```τ=V/F_in``` 
```K=-V/F_product``` 
```K_c=τ/(Kτ_c )```
```τ_i=τ```	(3)

The controllers were implemented in discrete form, using the following form for the ```k```-th value of the manipulated variable, ```u_k```, in terms of the proportional and integral gains, ```K_c``` and ```τ_i```, time step, ```Δt```, and error between setpoint and process value at timestep ```k```, ```e_k```.

```u_k=u_(k-1)+K_c (1+Δt/τ_i ) e_k```	(4)

### Tank network
The tank network is simulated using the function SimulateTankNetwork.m, which accepts the user specifications for the network, and outputs the measured values for the tank levels and all flow rates, for each tank, as well as the controller tuning and calculated steady state values for all parameters.

For the simple four tank network in the example code (mainFourTanks.m), the measured levels and flow rates are as shown in ```FourTankOutput.tif``` in the repo; note that 600 minutes of data were generated, for visual clarity, but 10 000 minutes of data were generated in the study described below.

The feed to the network enters tank 1, with product flow controlling level flowing to tank 2, and an underflow connection to tank 3. Tank 2’s product flow controls level, and leaves the system; the underflow feeds to tank 4. Tanks 3 and 4 do not have underflows; their product flows leave the system. More detail of the network can be found in mainFourTanks.m.

## References
Chen, D., Seborg, D.E., 2002. PI/PID controller design based on direct synthesis and disturbance rejection. Ind. Eng. Chem. Res. 41, 4807–4822. https://doi.org/10.1021/ie010756m

Farenzena, M., Trierweiler, J.O., 2009. LoopRank: A Novel Tool to Evaluate Loop Connectivity. IFAC Proc. Vol. 42, 964–969. https://doi.org/10.3182/20090712-4-TR-2008.00158

Lindner, B., Auret, L., Bauer, M., 2017a. Investigating the Impact of Perturbations in Chemical Processes on Data-Based Causality Analysis. Part 1: Defining Desired Performance of Causality Analysis Techniques. IFAC-PapersOnLine 50, 3269–3274. https://doi.org/10.1016/J.IFACOL.2017.08.463

Lindner, B., Auret, L., Bauer, M., 2017b. Investigating the Impact of Perturbations in Chemical Processes on Data-Based Causality Analysis. Part 2: Testing Granger Causality and Transfer Entropy. IFAC-PapersOnLine 50, 3275–3280. https://doi.org/10.1016/J.IFACOL.2017.08.620

Wakefield, B.J., Lindner, B.S., McCoy, J.T., Auret, L., 2018. Monitoring of a simulated milling circuit: Fault diagnosis and economic impact. Miner. Eng. 120, 132–151. https://doi.org/10.1016/j.mineng.2018.02.007

