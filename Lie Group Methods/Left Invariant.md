Multiply on left by group element and the vector field does not change

## Matrix differential equation
$\dot{X} = X \omega^\wedge$

## Lie algebra differential equation
$$\der{}{t}(\exp(\xi^\wedge) = X)$$
$$\der{}{t}(\exp(\xi^\wedge)) = \exp(\xi^\wedge)(\frac{I - \exp(-\ad{\xi^\wedge})}{\ad{\xi^\wedge}})\dot\xi^\wedge = \dot{X}$$
$$X^{-1}\dot{X} = (\frac{I - \exp(-\ad{\xi^\wedge})}{\ad{\xi^\wedge}})\dot\xi^\wedge = {\omega}^{\wedge}$$
$$\dot{\xi}^\wedge = \frac{\ad{\xi^\wedge}}{1 - \exp(-\ad{\xi^\wedge})}\omega^\wedge = J_l^{-1}(\xi)\omega^\wedge$$
Where $J_l$ is the left jacobian
$$\dot{\bar{X}} = \bar{X}\bar{\omega}$$
$$\eta = X^{-1}\bar{X}$$
$$\dot{\eta} = -X^{-1}\dot{X}X^{-1}\bar{X} + X^{-1}\dot{\bar{X}}$$
$$\dot\eta = -\omega^\wedge\eta + \eta\bar\omega^\wedge$$
$$\dot\eta = \eta\bar\omega^{\wedge} - \omega^{\wedge}\eta$$
Let $\omega = \bar\omega$
$$\eta^{-1}\dot{\eta} = \omega^\wedge - \eta^{-1}\omega^{\wedge}\eta$$
$$\eta^{-1}\dot\eta = \omega^{\wedge} - \Ad{\eta^{-1}}{\omega^\wedge}$$
$$\eta^{-1}\dot\eta = ((I - \Ad{\eta^{-1}}{})\omega)^\wedge$$
This means that $\eta^{-1}\dot\eta$ lives in the Lie algebra
Note: $\Ad{\eta - 1}{} = \exp(-\ad{\xi})$ and $\eta = \exp(\xi)$
$$\dot\xi = (\ad{\xi^{\wedge}}\omega)^\wedge$$
$$\dot{\xi} = \ad{\xi}\omega = -\ad{\omega}\xi$$
Therefore, the autonomous error dynamics are 
$$\dot{\xi} = -\ad{\omega}\xi$$
Which is LTI if $\omega$ is constant
BUT THIS IS ASSUMING $\omega = \bar\omega$  - This can be fixed with [[Log Linear Dynamic Inversion]]
