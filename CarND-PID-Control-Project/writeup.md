# CarND-Controls-PID
Self-Driving Car Engineer Nanodegree Program

---

## P, I, D components

P: the P gain is responsible for changing the steering wheel proportionally and negative to the error.
   So, for example, if the error of the position of the car are 3, the steering wheel will do a -P*3 movement. Because it always result in the same values for negative or positive errors, it tends to not stay stabilized on the set point line. Setting this gain very large are responsible for very sudden movements of the steering wheel.

D: derivative component, it is used to consider the change of the error. So, for example, when our P controller is running in a wave form, the derivative component can be used to keep our car in a movement closer to the set point curve (stable).

I: is responsible for the integral of the error of the position. It is used when our PD lead us trajectory far away from the goal, i.e, we have a bias between of our stable PD controller and our track.


## Tunning the hyperparameters

I started the tunning phase of the project by setting values to the P variable. As we saw on the classes, the P variable is responsible to get the output of our system around our set point, while in a zigzag form. In the beggining, I set the variable very large and then I got a very unstable wheel - because the proportional gain was very large, the steering values were basically with -1 or 1 input, i.e, the limits of the movement of the wheel.

So, I tried to decrease the P variable until I get a car that is capable to, at least, holds a track around the set point.

After that, I noticied that my car was not stabilizing itself in the proposed track. I tried to increase the gain of derivative - this was responsible to stopping the car from move in a wave form. I ended with a controller with P = 5, I=0, D=30.

Then, I noticied that the car was not able to handle some bends on the track and tried to increase the gain values for P and D. One thing interesting here is that my simulator was running locally, while my PID controller was running remotely 100km far away. I did realize that the delay between the controller and the car was difficulting the tunning of the hyperparameters a lot!

Once running it locally, I increased the values and got the car driving better. I tried to  put I=0.4 but I noticied that the variable created a lot of "momentum" in the car error, doing the car go far from the lane boundaries. I then tried I=0.1 and noticied that it was better, but the problem was still there.

Finally, I ended my parameters with P=7, I=0, D=40. Yet, its necessary to remove some unstability around the set point position. I tried to increase the Derivative gain but it seems that it can very dangerous (specially for a human inside a car) due to very strong movements that the car in some moments.

Twiddle can be used to try to fine tunning the PID controller as an improvement of this project.
