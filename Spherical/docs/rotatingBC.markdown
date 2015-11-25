Spherical Underworld
====================

Making a spherical mesh
---------------------


Solving on a spherical mesh
---------------------------
We seek to solve Stokes Flow on the Spherical mesh without rewriting the equations in new coordinate (i.e. r, theta, phi )
Engelman et al. (1982) describe a technique for implementing boundary condition which no longer align with the x-y-z axis. We use this method within the existing code to solve on the Spherical mesh with the equation still in x-y-z

Consider this 2D example of how the Stokes momentum equation works for a single boundary node.

$K\boldsymbol{u}=\boldsymbol{f}$

where $\boldsymbol{u}$ is specified in x-y. Let's take this node to be on a boundary surface that doesn't align with the x or y axis and we wish for the node to have a constraint velocity for it normal component and uncontraint for the tangential component i.e. 

$u_{n}=\boldsymbol{n}\cdot\boldsymbol{u} = n_{x}u_{x} + n_{y}u_{y} = \phi$

More generally we can rotate a given $(u_{x}, u_{y})$ into $(u_{n}, u_{t})$ via 

$\boldsymbol{ u^{'}}
=
\begin{bmatrix}
u_{n}\\
u_{t} 
\end{bmatrix}
= 
\begin{bmatrix}
n_{x} & n_{y}\\ 
t_{x} & t_{y}
\end{bmatrix}
\begin{bmatrix}
u_{x}\\
u_{y}
\end{bmatrix}
=
R^{T}\boldsymbol{u}$

Where $R$ is a rotation matrix such that $R^{T}=R^{-1}$

Therefore we can rewrite the l.h.s of the momentum equation with $\boldsymbol{u'}$

$K\boldsymbol{u}=KRR^{T}\boldsymbol{u}=KR\boldsymbol{u'}$

This method translates normal and tangent velocity BC into the x-y axis-aligned momentum equation. Upon these rotations incorrect $\boldsymbol{u'}$ can be applied to x-y axis-aligned momentum equations. To fix this problem we can rewrite the momentum equation in its normal and tangential component, i.e. rotate the momentum equation to be 

$M^{'}
=
\begin{bmatrix}
M_{n}\\
M_{t} 
\end{bmatrix}
= 
\begin{bmatrix}
n_{x} & n_{y}\\ 
t_{x} & t_{y}
\end{bmatrix}
\begin{bmatrix}
M_{x}\\
M_{y}
\end{bmatrix}
=
R^{T}M$

Applying these methods of rotating the momentum equation by PRE-MULTIPLYING by $R^{T}$ and rewriting PRE-MULTIPLYING $\boldsymbol{u'}$ by $R$:

$\begin{bmatrix}
R^{T}K & R^{T}G\\ 
G^{T} & 0
\end{bmatrix}
\begin{bmatrix}
R\boldsymbol{u'}\\
p
\end{bmatrix}
=
\begin{bmatrix}
R^{T}\boldsymbol{f}\\
0
\end{bmatrix}$

When can rewrite this so the rotation are only applied to the assembled stiffness matricies and assembled r.h.s vectors. This minimises the changes to the existing code structure because only nodal $R$ need to be calculated for nodes that have dofs not defined in the x-y-z axis. In Underworld's implementation elemental $R$ are calculated for all nodes belonging to a special set and applied as follows

$\begin{bmatrix}
R^{T}KR & R^{T}G\\ 
G^{T}R & 0
\end{bmatrix}
\begin{bmatrix}
\boldsymbol{u'}\\
p
\end{bmatrix}
=
\begin{bmatrix}
R^{T}\boldsymbol{f}\\
0
\end{bmatrix}$


References
----------
* Engelman et al. THE IMPLEMENTATION OF NORMAL AND/OR TANGENTIAL BOUNDARY CONDITIONS IN FINITE ELEMENT CODES FOR INCOMPRESSIBLE FLUID FLOW, Int. JOURNAL FOR NUMERICAL METHODS IN FLUIDS, VOL. 2, 225-238 (1982)
