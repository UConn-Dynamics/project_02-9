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
