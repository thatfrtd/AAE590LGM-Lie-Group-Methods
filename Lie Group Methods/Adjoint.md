
Adjoint is the linear map mapping from Lie algebra to tangent space at a group element
$\Ad{X}{\cdot}:T_I G \to T_X G$

$v^\wedge = V \in \mathfrak{g}$
$$\Ad{X}{v} = Xv^\wedge X^{-1} = (v')^{\wedge}$$
Want linear operator for adjoint $v' = \Ad{X}{v}$

## Special Orthogonal Group
Note: for $SO(2)$, $SO(3)$ the adjoint is rotation matrix itself

Recover ____ transformation

$$\omega_b = R^{bw}\omega_w = \Ad{X}{\omega_w}$$
$$\omega_b^{\wedge} = R\omega_w^{\wedge} R^T = \Ad{R}{\omega_w}^\wedge$$
## $SE(2)$
$$X = \Mtwo{R}{p}{}{1}, \: X^{-1} = \Mtwo{R^T}{-R^Tp}{}{1}$$
$$\Ad{X}{v} = Xv^\wedge X^{-1} = (v')^{\wedge} =  \Mtwo{R}{p}{}{1}\Mtwo{\omega^\wedge}{u}{}{}\Mtwo{R^T}{-R^Tp}{}{1} = \Mtwo{R}{p}{}{1}\Mtwo{\omega^\wedge R^T}{-\omega^\wedge R^Tp + u}{}{} = \Mtwo{\omega^\wedge}{-\omega^\wedge p + Ru}{}{} = \Mtwo{1}{-p}{}{}\Vthree{}{}{}$$
$\Ad{X}{} = \Mtwo{R}{p^\wedge}{}{1}$

