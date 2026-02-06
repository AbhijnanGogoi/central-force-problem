# central-force-problem

Documentation preparation in progress.

To run demo code:
`cd two-body-problem`
`sh demo.sh`

Output will be in `two-body-problem/demo` directory


## Description of the algorithm:
The classical two body central force problem is one of those problems in classical mechanics that have been solved analytically. However, the three-body and general n-body problems are not yet solved analytically. One of the ways of obtaining the particular solutions to such problems is the use of computer simulations.

The code for this simulation of the classical two-body central force problem (which I aim to extend to general n-body systems) uses the Hamiltonian formalism and tries to obtain the time curve of the system in the phase space. However, integrating with respect to time comes with its problems (which may have already been addressed in existing algorithms). One of these include the uneven step sizes with each iteration in the position space due to the changing velocities of the objects in the force field. These objects can move with trajectories that do not show much resemblance with a straight line. In such cases, the time step $\Delta t$ must be decreased further so as to ensure that the trajectory obtained by the simulation resembles with the actual trajectory followed by the bodies. However, this increases the computational load if a large part of the actual trajectory does not show much curvature and can be integrated with a larger time step to give an acceptable result. Hence, this simulation keeps the change in position $|\Delta \vec r|$ constant in each iteration. This will ensure that the trajectory is captured by our simulation close to the actual trajectory of the objects.

Although keeping track of time is not necessary as the velocity of the particles can be obtained from the kinetic energy of the system (which can be calculated from the total energy and the potential energy which can be obtained from the central force law), the simulation code is updating time along with the coordinates and canonical momenta of the system.

### The Hamilton's equations of motion:
The coordinates chosen for the simulation are the distance between the two bodies (given as $r$) and the orientation of the position vector of object 2 with respect to object 1 (given as $\theta$).

The canonical momenta of these coordinates are as follows:
$$
r \to p_r \\
\theta \to p_\theta
$$

The Hamiltonian of the system obtained is as follows:
$$H = \frac{p_r^2}{2\mu} + \frac{p_\theta^2}{2\mu r^2} + V(r)$$
where $V(r)$ is the potential energy function of the system and $\mu$ is the reduced mass of the system given by:
$$
\mu = \frac{m_1 m_2}{m_1 + m_2}
$$

In our simulation, we are taking $V(r) = kr^n$ where $n$ is any integer.

So, the Hamiltonian equations of motion:
$$
\dot\theta = \frac{p_\theta}{\mu r^2} \\
\dot r = \frac{p_r}{\mu} \\
\dot p_\theta = 0\\
\dot p_r = \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr}
$$


### The equations for updating the coordinates and momenta with respect to $|\Delta \vec r|$

We define a differential $dx^2 = dr^2 + r^2 d\theta^2$.
This gives:
$$
\left(\frac{dx}{dt}\right)^2 = \left(\frac{dr}{dt}\right)^2 + r^2 \left(\frac{d\theta}{dt}\right)^2
$$

Solving for $dx/dt$ gives:
$$
\frac{dx}{dt} = \frac{1}{\mu r} \sqrt{r^2 p_r^2 + p_\theta^2}
$$

Since we want to iterate with respect to change in $x$ and not $t$, we parameterise the coordinates, momenta and time using $x$ whose second-order Taylor approximation about $x_0$ (which is the current value of $x$) is as follows:
$$
r \equiv r(x) = r(x_0) + (x-x_0) r'(x_0) + \frac{(x-x_0)^2}{2} r''(x_0)
$$

$$
\theta \equiv \theta(x) = \theta(x_0) + (x-x_0) \theta'(x_0) + \frac{(x-x_0)^2}{2} \theta''(x_0)
$$

$$
p_r \equiv p_r(x) = p_r(x_0) + (x-x_0) p_r'(x_0) + \frac{(x-x_0)^2}{2} p_r''(x_0)
$$

$$
p_\theta \equiv p_\theta(x) = p_\theta(x_0) + (x-x_0) p_\theta'(x_0) + \frac{(x-x_0)^2}{2} p_\theta''(x_0)
$$

$$
t \equiv t(x) = t(x_0) + (x-x_0) t'(x_0) + \frac{(x-x_0)^2}{2} t''(x_0)
$$

where $f'$ denotes the derivative of $f$ with respect to $x$. A second-order and not a first-order Taylor approximation is used as the value of the change in the parameter (e.g. $\Delta r$) vanishes to zero for $p_r = 0$ and $p_\theta = 0$ simultaneously even when a non-zero force is acting on the objects.

The change in the parameters in a single iteration is obtained as:
$$
\Delta t = t'(x_0) \Delta x + \frac{\Delta x^2}{2} t''(x_0)
$$

$$
\Delta r = r'(x_0) \Delta x + \frac{\Delta x^2}{2} r''(x_0)
$$

$$
\Delta \theta = \theta'(x_0) \Delta x + \frac{\Delta x^2}{2} \theta''(x_0)
$$

$$
\Delta p_r = p_r'(x_0) \Delta x + \frac{\Delta x^2}{2} p_r''(x_0)
$$

$$
\Delta p_\theta = 0
$$

where $\Delta x$ is equal to $x-x_0$ which is kept constant throughout the simulation. $\Delta p_\theta$ is zero as $p_\theta$ is a conserved quantity as seen from the Hamilton's equations of motion.

So in each iteration, the program calculates the first and second-order derivatives $r'$, $r''$, $\theta'$, $\theta''$, $p_r'$, $p_r''$, $t'$ and $t''$ at $x=x_0$ (that is the current value of $x$) using which it obtains the change in the parameters required and the values of the time, coordinates and the momenta in the next iteration.


### The first and second-order derivatives of the parameters with respect to x
#### Case 1: (atleast one of $p_r$ and $p_\theta$ is non-zero)
Using $dx/dt$ from the above, we know that:
$$
t'(x_0) = \frac{\mu r}{\sqrt{r^2 p_r^2 + p_\theta^2}}
$$

Hence, the first-order derivatives with respect to $x$ are obtained from the Hamilton's equations using the chain rule with $dt/dx$ as:

$$
r'(x_0) = \frac{r p_r}{\sqrt{r^2 p_r^2 + p_\theta^2}}
$$

$$
\theta'(x_0) = \frac{p_\theta}{r\sqrt{r^2 p_r^2 + p_\theta^2}}
$$

$$
p_r'(x_0) = \frac{\mu r}{\sqrt{r^2 p_r^2 + p_\theta^2}} \left[ \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr} \right]
$$


Again taking the derivative with respect to $x$ and substituting the above results gives the following second-order derivatives:

$$
r''(x_0) = rp_r \left[ \frac{p_r}{\lambda} + \frac{\mu r}{p_r \lambda}\left( \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr} \right) - \frac{\mu r^3 p_r}{\lambda^2} - \frac{r^2 p_r^3}{\lambda^2} \right]
$$

$$
\theta''(x_0) = \frac{-p_\theta}{r} \left[ \frac{p_r}{\lambda} + \frac{r^2 p_r^3}{\lambda^2} + \frac{\mu r^3 p_r}{\lambda^2}\left( \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr} \right) \right]
$$

$$
p_r''(x_0) = \mu r \left( \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr} \right) \left[ \frac{p_r}{\lambda} - \frac{r^2 p_r^3}{\lambda^2} - \frac{\mu r^3 p_r}{\lambda^2} \left( \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr} \right) - \frac{3p_\theta^2 p_r}{\mu r^3 \lambda} - (n-1)\frac{p_r}{\lambda} \frac{dV}{dr} \right]
$$

$$
t''(x_0) = \mu r \left[ \frac{p_r}{\lambda} - \frac{r^2 p_r^3}{\lambda^2} - \frac{p_r p_\theta^2}{\lambda^2} + \frac{\mu r^3 p_r}{\lambda^2}\frac{dV}{dr} \right]
$$

where $\lambda = r^2 p_r^2 + p_\theta^2$, and $V(r) = k r^n$.

#### Case 2: (both $p_r$ and $p_\theta$ are zero)

In this case, it suggests that in the current state of the system (at $x=x_0$) the system has zero velocity. Since $p_\theta$ is a conserved quantity, we get $\Delta \theta = 0$. Hence, $|\Delta r| = |\Delta x|$. Since the direction of change of $r$ depends on the direction of the central force ($-dV/dr$), we get:

if $-dV/dr > 0$, $\Delta r = \Delta x$
if $-dV/dr = 0$, $\Delta r = 0$
if $-dV/dr < 0$, $\Delta r = -\Delta x$

Use the equations of motion for constant acceleration gives us:
$$\Delta r = \frac{dr}{dt} \Delta t + \frac{1}{2} \frac{d^2r}{dt^2} \Delta t^2$$
which gives us (for $-dV/dr \neq 0$):

$$
\Delta t = \sqrt{\frac{2 \mu \Delta x}{|-dV/dr|}}
$$

For $-dV/dr = 0$, we see that the system is in equilibrium and stationary, so there is no change in any parameters which makes us free to choose $\Delta t$.

From Hamilton's equations of motion, we know that:
$$\dot p_r = \frac{p_\theta^2}{\mu r^3} - \frac{dV}{dr}$$
which gives us:
$$\Delta p_r = \frac{-dV}{dr} \Delta t$$

