$$X^{-1}\dot{X} = \omega^\wedge \Leftrightarrow \dot{X} = X\omega^\wedge$$
$$\dot{\xi} = -\text{ad}_\bar{\omega}\xi + J_{l(\xi)}^{-1}\tilde\omega$$
If $\tilde\omega = J_{l(\xi)} BK\xi$ so $\tilde{\omega}$ is basically the control action
then $\dot{\xi} = A\bar\omega\xi + BK\xi = (A\bar\omega + BK)\xi$ which is nice and linear
## Underactuated Control
If you have a constraint that your control $\tilde\omega$ can only act in certain directions then even if you design the feedback matrix $K$ to satisfy this, after the transformation by $J_{l(\xi)}$ this is not likely to be satisfied. One solution is to bound the term and call it a disturbance. If you bound the term and analyze your system properly you can show that the controller still works well.

## Paper
[[Log-Linear Dynamic Inversion Control With Provable Safety Guarant.pdf]]