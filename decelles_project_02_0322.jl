### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ f17103ea-06bf-11f1-a2b0-79e68ed152eb
md"""# Project_02 - Multibody kinematic modeling

![Dual slider kinematics project](https://raw.githubusercontent.com/cooperrc/me5180-project_02/refs/heads/main/dual-slider.svg)

In this project, a rigid bar is connected to two sliding pistons along
the diagonal tracks. As the pistons move along the tracks, the rigid bar rotates at a constant rate, $\dot{\theta}_3 = 2~rad/s$. The figure above has three _relative_ ccoordinate systems that move with the bodies:

1. $x_1-y_1-$ describes piston 1 position and orientation, $\theta_1$
2. $x_2-y_2-$ describes piston 2 position and orientation, $\theta_2$
3. $x_3-y_3-$ describes the rigid bar position and orientation, $\theta_3$

Each of the pistons are on tracks at $\pm 45^o$ and the rotating rigid
bar is 10 cm. The hinges are mounted to the center of the pistons
connecting the ends of the rigid bar. 
 
In this project, you need to 

1. determine constraint equations $C(\mathbf{q},~t)$
2. solve for the velocities, $\dot{q}$ and accelerations, $\ddot{q}$
3. visualize the motion of the system as the rigid bar goes through at least one full rotation
"""

# ╔═╡ 0d9be664-d7c5-4084-add2-25e5418742d6
md"""# 1. Determine the Constraint equations C(q,t)

Knowns:

1. $x_1-y_1-$ describes piston 1 position and orientation, $\theta_1$
2. $x_2-y_2-$ describes piston 2 position and orientation, $\theta_2$
3. $x_3-y_3-$ describes the rigid bar position and orientation, $\theta_3$

$\dot{\theta_3}= 2~rad/s$

$q=[x_1, x_2, \theta_3]$

$L=10~cm=0.1~m$


Setting up Constraints:

Constraint 1

$\mathbf{r}(t)_{1} = \mathbf{r}'+{x}_{1}i'+ {y}_{1}j'$

$\mathbf{r}(t)_{1} = \mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'$

$\mathbf{r}(t)_{1} = \mathbf{r}'^G +{m}_{x}i'^0 + {m}_{y}j'^0$

$\mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'= \mathbf{r}'^G +{m}_{x}i'^0 + {m}_{y}j'^0$

$\mathbf{q}(t)_{1} = 
\begin{bmatrix}
\mathbf{r} \\
{\theta}'
\end{bmatrix}=\begin{bmatrix}
{r}_{x}' \\
{r}_{y}' \\
{\theta}'
\end{bmatrix}$

$C(\mathbf{q},~t)=\mathbf{0}$

$\mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'- \mathbf{r}'^G - {m}_{x}i'^0 - {m}_{y}j'^0=\mathbf{0}$

${r}_{x}'-(\frac{w}{2}i'+\frac{t}{2}j')i^0-0-{m}_{x}=0$
${r}_{y}'-(\frac{w}{2}i'+\frac{t}{2}j')j^0-0-{m}_{y}=0$

$i'=cos{\theta}'i^0+sin{\theta}'j^0$
$j'=-sin{\theta}'i^0+cos{\theta}'j^0$

$\begin{bmatrix}
i' \\
j'
\end{bmatrix}=\begin{bmatrix}
cos{\theta}'~sin{\theta}' \\
-sin{\theta}'~cos{\theta}'
\end{bmatrix}\begin{bmatrix}
i^0 \\
j^0
\end{bmatrix}$

$\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{1}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{1}=\mathbf{0}$

$\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{1}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{1}=\mathbf{0}$

$\mathbf{C}(\mathbf{q},~t)=\mathbf{0}=\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{1}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{1}$

$\frac{d}{dt}C(\mathbf{q},~t)=velocity~constraints$

$\frac{d}{dt}(\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{2}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{2})$

$\mathbf{v}'+\frac{d}{dt}(\mathbf{A}')\mathbf{r}'_{1}+\mathbf{A}'\frac{d}{dt}(\mathbf{r}'_{1})+0+0$

$\frac{d{\theta}}{dt}\frac{d}{d{\theta}}(\begin{bmatrix}
cos{\theta}'~-sin{\theta}' \\
sin{\theta}'~cos{\theta}'
\end{bmatrix})$

$\mathbf{0}=\mathbf{v}'+ \dot{\theta}'\begin{bmatrix}
-sin{\theta}'~-cos{\theta}' \\
cos{\theta}'~-sin{\theta}'
\end{bmatrix}\mathbf{r}'_{1}-\mathbf{v}^G-\mathbf{\dot{\theta}}^0\begin{bmatrix}
-sin{\theta}'~-cos{\theta}' \\
cos{\theta}'~-sin{\theta}'
\end{bmatrix}\mathbf{r}'_{1}$

$\frac{d}{dt}(C(\mathbf{q},~t))=\mathbf{0}$

$\frac{d\mathbf{C}}{d\mathbf{q}}\mathbf{\dot{q}}+\frac{d\mathbf{C}}{dt}=0$

$\mathbf{r}'+\mathbf{A}'\begin{bmatrix}
\frac{-w}{2} \\
\frac{t}{2}
\end{bmatrix}=0=C(\mathbf{q},~t)$

$\mathbf{r}'_{x}+\frac{-w}{2}cos{\theta}'-\frac{t}{2}sin{\theta}'=0$

$\mathbf{r}'_{y}+\frac{-w}{2}sin{\theta}'+\frac{t}{2}cos{\theta}'=0$

$\frac{dC'}{d{r}'_{x}}=1$

$\frac{dC'}{d{r}'_{y}}=0$

$\frac{dC'}{d{\theta}^T}=\frac{w}{2}sin{\theta}'-\frac{t}{2}cos{\theta}'$

$\frac{dC^2}{d{r}'_{x}}=0$

$\frac{dC^2}{d{r}'_{y}}=1$

$\frac{dC'}{d{\theta}'}=\frac{-w}{2}cos{\theta}'-\frac{t}{2}sin{\theta}'$

$\frac{dC'}{dq'}=\begin{bmatrix}
1~0~\mathbf{A}{\theta}'\frac{-w}{2} \\
0~1~\mathbf{A}{\theta}'\frac{t}{2} 
\end{bmatrix}$

$\frac{dC'}{dq'}\begin{bmatrix}
{v}'_{x} \\
{v}'_{y} \\
{\dot{\theta}}'
\end{bmatrix}=\begin{bmatrix}
\frac{dC'_1}{dt} \\
\frac{dC'_2}{dt}
\end{bmatrix}=0$

${C_3}={r_y}'-2t=0$

$C(\mathbf{q},~t)=\begin{bmatrix}
\mathbf{r}'+\mathbf{A}\begin{bmatrix}
\frac{-w}{2} \\
\frac{t}{2}
\end{bmatrix} \\
\mathbf{r_y}'-2t 
\end{bmatrix}=\mathbf{0}$

$\frac{d\mathbf{C'}}{d\mathbf{q'}}=\begin{bmatrix}
1~0~\mathbf{A}{\theta}'\frac{-w}{2} \\
0~1~\mathbf{A}{\theta}'\frac{t}{2} \\
0~1~0
\end{bmatrix}\begin{bmatrix}
{v}'_{x} \\
{v}'_{y} \\
{\dot{\theta}}'
\end{bmatrix}=\begin{bmatrix}
0 \\
0 \\
2
\end{bmatrix}$

$\mathbf{C}(\mathbf{q}_{\mathbf{q}},\mathbf{q})=-\mathbf{C_t}$

$\frac{d\mathbf{C}}{d\mathbf{q}}=\mathbf{C_q}$

$\frac{d\mathbf{C}}{dt}=\mathbf{C}_{t}$

$\frac{d^2\mathbf{C}}{dt^2}=\mathbf{C}_{tt}$


Constraint 2


$\mathbf{r}(t)_{2} = \mathbf{r}'+{x}_{2}i'+ {y}_{2}j'$

$\mathbf{r}(t)_{2} = \mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'$

$\mathbf{r}(t)_{2} = \mathbf{r}'^G +{m}_{x}i'^0 + {m}_{y}j'^0$

$\mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'= \mathbf{r}'^G +{m}_{x}i'^0 + {m}_{y}j'^0$

$\mathbf{q}(t)_{1} = 
\begin{bmatrix}
\mathbf{r} \\
{\theta}'
\end{bmatrix}=\begin{bmatrix}
{r}_{x}' \\
{r}_{y}' \\
{\theta}'
\end{bmatrix}$

$C(\mathbf{q},~t)=\mathbf{0}$

$\mathbf{r}'-\frac{w}{2}i'+\frac{t}{2}j'- \mathbf{r}'^G - {m}_{x}i'^0 - {m}_{y}j'^0=\mathbf{0}$

${r}_{x}'-(\frac{w}{2}i'+\frac{t}{2}j')i^0-0-{m}_{x}=0$
${r}_{y}'-(\frac{w}{2}i'+\frac{t}{2}j')j^0-0-{m}_{y}=0$

$i'=cos{\theta}'i^0+sin{\theta}'j^0$
$j'=-sin{\theta}'i^0+cos{\theta}'j^0$

$\begin{bmatrix}
i' \\
j'
\end{bmatrix}=\begin{bmatrix}
cos{\theta}'~sin{\theta}' \\
-sin{\theta}'~cos{\theta}'
\end{bmatrix}\begin{bmatrix}
i^0 \\
j^0
\end{bmatrix}$

$\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{2}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{2}=\mathbf{0}$

$\mathbf{C}(\mathbf{q},~t)=\mathbf{0}=\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{2}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{2}$

$\frac{d}{dt}C(\mathbf{q},~t)=velocity~constraints$

$\frac{d}{dt}(\mathbf{r}'+\mathbf{A}'\mathbf{r}'_{2}-\mathbf{r}^G-\mathbf{A}^0\mathbf{r}^0_{2})$

$\mathbf{v}'+\frac{d}{dt}(\mathbf{A}')\mathbf{r}'_{2}+\mathbf{A}'\frac{d}{dt}(\mathbf{r}'_{2})+0+0$

$\frac{d{\theta}}{dt}\frac{d}{d{\theta}}(\begin{bmatrix}
cos{\theta}'~-sin{\theta}' \\
sin{\theta}'~cos{\theta}'
\end{bmatrix})$

$\mathbf{0}=\mathbf{v}'+ \dot{\theta}'\begin{bmatrix}
-sin{\theta}'~-cos{\theta}' \\
cos{\theta}'~-sin{\theta}'
\end{bmatrix}\mathbf{r}'_{2}-\mathbf{v}^G-\mathbf{\dot{\theta}}^0\begin{bmatrix}
-sin{\theta}'~-cos{\theta}' \\
cos{\theta}'~-sin{\theta}'
\end{bmatrix}\mathbf{r}'_{2}$

$\frac{d}{dt}(C(\mathbf{q},~t))=\mathbf{0}$

$\frac{d\mathbf{C}}{d\mathbf{q}}\mathbf{\dot{q}}+\frac{d\mathbf{C}}{dt}=0$

$\mathbf{r}'+\mathbf{A}'\begin{bmatrix}
\frac{-w}{2} \\
\frac{t}{2}
\end{bmatrix}=0=C(\mathbf{q},~t)$

$\mathbf{r}'_{x}+\frac{-w}{2}cos{\theta}'-\frac{t}{2}sin{\theta}'=0$

$\mathbf{r}'_{y}+\frac{-w}{2}sin{\theta}'+\frac{t}{2}cos{\theta}'=0$

$\frac{dC'}{d{r}'_{x}}=1$

$\frac{dC'}{d{r}'_{y}}=0$

$\frac{dC'}{d{\theta}^T}=\frac{w}{2}sin{\theta}'-\frac{t}{2}cos{\theta}'$

$\frac{dC^2}{d{r}'_{x}}=0$

$\frac{dC^2}{d{r}'_{y}}=1$

$\frac{dC'}{d{\theta}'}=\frac{-w}{2}cos{\theta}'-\frac{t}{2}sin{\theta}'$

$\frac{dC'}{dq'}=\begin{bmatrix}
1~0~\mathbf{A}{\theta}'\frac{-w}{2} \\
0~1~\mathbf{A}{\theta}'\frac{t}{2} 
\end{bmatrix}$

$\frac{dC'}{dq'}\begin{bmatrix}
{v}'_{x} \\
{v}'_{y} \\
{\dot{\theta}}'
\end{bmatrix}=\begin{bmatrix}
\frac{dC'_1}{dt} \\
\frac{dC'_2}{dt}
\end{bmatrix}=0$

${C_3}={r_y}'-2t=0$

$C(\mathbf{q},~t)=\begin{bmatrix}
\mathbf{r}'+\mathbf{A}\begin{bmatrix}
\frac{-w}{2} \\
\frac{t}{2}
\end{bmatrix} \\
\mathbf{r_y}'-2t 
\end{bmatrix}=\mathbf{0}$

$\frac{d\mathbf{C'}}{d\mathbf{q'}}=\begin{bmatrix}
1~0~\mathbf{A}{\theta}'\frac{-w}{2} \\
0~1~\mathbf{A}{\theta}'\frac{t}{2} \\
0~1~0
\end{bmatrix}\begin{bmatrix}
{v}'_{x} \\
{v}'_{y} \\
{\dot{\theta}}'
\end{bmatrix}=\begin{bmatrix}
0 \\
0 \\
2
\end{bmatrix}$

$\mathbf{C}(\mathbf{q}_{\mathbf{q}},\mathbf{q})=-\mathbf{C_t}$

$\frac{d\mathbf{C}}{d\mathbf{q}}=\mathbf{C_q}$

$\frac{d\mathbf{C}}{dt}=\mathbf{C}_{t}$

$\frac{d^2\mathbf{C}}{dt^2}=\mathbf{C}_{tt}$
"""

# ╔═╡ 825b962e-51b8-424a-8378-8fe4675e1033
md"""# 2. Solve for Velocities $\dot{q}$ and accelerations $\ddot{q}$

"""

# ╔═╡ eed36891-6d9c-4a69-9285-0bc6df359dcd
md"""# 3. Visualize the Motion of the System as a Rigid Bar Goes through at least one Full Rotation

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╟─f17103ea-06bf-11f1-a2b0-79e68ed152eb
# ╠═0d9be664-d7c5-4084-add2-25e5418742d6
# ╠═825b962e-51b8-424a-8378-8fe4675e1033
# ╠═eed36891-6d9c-4a69-9285-0bc6df359dcd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
