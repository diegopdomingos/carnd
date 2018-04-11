
##The model

The model used on MPC can be splitted on two parts: the cost due high actuators, changes, etc and the constraints due the phisical model. 

#Total Cost

The total cost of the MPC is calculated and stored on the fg[0] vector position. In my model I considered the following:

For errors in general:
 
* Cross track error (cte(t))
* Car Direction error (epsi(t))
* Car velocity error (v-ref_v(t))

For actuators:

* Steering actuator (delta)
* Accelerator actuator (a)
* Change in steering actuator (delta(t+1)-delta(t))
* Change in accelerator actuator (a(t+1)-a(t)

For steering actuator, accelerator actuator and change in steering actuator there is a weight to be determined, defined as delta_weight, a_weight and delta_var_weight float variables, respectively.

#Update equations

The model was implemented in the same way we learned on Udacity Self-Driving Car classes, so:

x1   = (x0 + v0*cos(psi0)*dt);
y1   = (y0 + v0*sin(psi0)*dt);
psi1 = (psi0 - v0*delta0/Lf*dt);
v1   = (v0 + a0*dt);
cte1 = ((f0 - y0) + (v0*sin(epsi0)*dt));
epsi1= ((psi0 - psides0) - v0*delta0/Lf*dt);


#Timestep length and time

The timestep length chosen was N=10 and the time for each timestep was 0.1. So, for each timestep, we have the same amount of delay time of actuators. This will be important since this is required and cannot be changed due the way that I managed to consider the actuator delay, explained further in details in this document.

The N and dt values were chosen considering:

- 0.1s is the same ammount of actuators delay. So we can consider the our actuators will work in every second timestep of our model (t+1)
- each interation we can determine the cost function based in 1 second of the movement, which is reasonable for a vehicle running in about 20 m/s. So, in each interaction we are dealing with 20 meters forward.
- the number of calculations to avoid a lot of processing and delay.
- the scale of roads (we are dealing with 20 meters forward)

##Polynomial Fitting and MPC Preprocessing

# Process waypoints

The waypoints are transformed to the car origin using the defined function transformWaypointsToCar(). This functions basically creates two rows: the row 0 defines the transformed x points (ptsx) and the row 1 defines the transformed y points (ptsy), in relation of the car position (px,py).

# Polynomial fit

After the transformation, we now need to fit our points (transformed waypoints) on a 3rd degree curve. This is done by the polyfit function.

# State

Since our x0, y0 and psi0 are zero, i.e, we consider the start of curve in the origin and the curve is in the form A + Bx + Cx^2 + Dx^3, we have the following state:

- x0   = 0
- y0   = 0
- psi0 = 0
- cte0 = A
- epsi0 = -atan(B);

# Model Predictive Control with Latency

Because the time for each step in the model is dt=0.1s, we have one timestep of delay (100ms). Because of that, we can consider always the actuator of one step forward, or, how it was implemented:

x1.push_back(solution.x[delta_start+1]);
x1.push_back(solution.x[a_start+1]);

In this case, we are considering always the delta_start of the one timestep forward, which is the actuators for the state 100ms forward.

