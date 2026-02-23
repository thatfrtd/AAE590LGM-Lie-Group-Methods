
$v$ = wind speed
$u$ = airspeed

Dynamics
$$\begin{gather*} 
\dot{x} = u\cos(\theta) + v \\
\dot{y} = u\sin(\theta) \\
\dot{\theta} = \omega
\end{gather*}$$
Group element
$X = \Mtwo{R(\theta)}{\vec{p}(x,y)}{}{1} = \Mthree{c\theta}{-s\theta}{x}{s\theta}{c\theta}{y}{0}{0}{1}$

$$\dot{X} = \Mthree{-\omega s\theta}{-\omega s\theta}{u c\theta + v}{\omega c\theta}{-\omega s\theta}{u s \theta}{0}{0}{0} = V_r X + XV_l$$
$V_r$ is [[right invariant]] vector field (world)
$V_l$ is [[left invariant]] vector field (body)
System is [[mixed invariant]]

The airspeed lives in the tangent space at the current state $X$.

$V_r = \Mthree{0}{0}{v}{0}{0}{0}{0}{0}{0}$
$V_l = \Mthree{0}{-\omega}{u}{\omega}{0}{0}{0}{0}{0}$


