# **Power Flow**


## Introduction

The power flow is the analysis of the flow of electric power, voltages, and currents in normal steady-state operation. It calculates the voltage magnitude and phase angle of each bus from which real and reactive flow on all lines as well as the associated equipment losses can be computed.
The mathematical model of the power flow problem is formulated as a set of nonlinear equations. Due to the nonlinear nature of the problem formulation, power flow does not have an explicit solution, but rather a numerical iterative approaches are used to obtain the solution. Input data to a power flow problem consist of bus data, transmission line data, and transformer data.

## Network Equations

The relation between all the bus current injections and bus voltage is given by the following node equation: 
```math
I_{i}=\sum_{j=1}^{n} Y_{ij}V_{j} \;\;\;\;\;\; i=1,2,...,n
```
where  $`n`$ is the number of buses in the network. $`I_{i}`$ and $`V_{j}`$ are injected current at bus $`i`$ and voltage at bus $`j`$. $`Y_{ij}`$ are the elements of the admittance matrix **Y**. Diagonal elements $`Y_{ii}`$ are equal to the sum of all admittances of all devices incident to the bus $`i`$. Off-diagonal elements $`Y_{ij}`$ are equal to the **negative** of the sum of the admittances that are joining buses $`i`$ and $`j`$. In case that there is shift transformer at the bus, $`Y_{ij}`$ should be calculated as explained in the branch section.

In the power system, complex voltage and current values are unknown, but rather real power injections at the generator buses and voltage magnitude setpoint as well as complex power (S) consumed by the load.
The relation between injected current and power at the node is given as:
```math
S_{i}=P_{i}+jQ_{i}=V_{i}I^*_{i}
```
Mode information on the sign convention used can be found [here](/ComponentLib/Bus/README.md).

Now the complex power flow equations are given as:

```math
P_{i}+jQ_{i}=V_{i}\sum_{j=1}^{n} Y_{ij}^*V_{j}^* \;\;\;\;\;\; i=1,2,...,n
```

Considering:
```math
V_{i}=\vert V_{i} \vert e^{j\theta_{i}}
```

```math
Y_{ij}=\vert Y_{ij} \vert e^{j\psi_{ij}}=G_{ij}+jB_{ij}
```

```math
\theta_{ij}=\theta_{i}-\theta_{j}
```

The complex set can be rewritten as two sets of real power balance equations:

```math
P_{i}= \vert V_{i} \vert \sum_{j=1}^{n}\vert Y_{ij} \vert \vert V_{j} \vert \cos (\theta_{ij}-\psi_{ij})
```

```math
Q_{i}= \vert V_{i} \vert \sum_{j=1}^{n} \vert Y_{ij} \vert \vert V_{j} \vert \sin (\theta_{ij}-\psi_{ij})
```

or

```math
P_{i}= \vert V_{i} \vert \sum_{j=1}^{n}\vert V_{j} \vert (G_{ij}\cos\theta_{ij}+B_{ij}\sin\theta_{ij})
```

```math
Q_{i}= \vert V_{i} \vert \sum_{j=1}^{n} \vert V_{j} \vert (G_{ij}\sin\theta_{ij}-B_{ij}\cos\theta_{ij})
```


## Example Power Flow Calculaion






<div align="center">
   <img align="center" src="https://gitlab.pnnl.gov/gridkit/gridkit/-/raw/my_develop_new/Documentation/Figures/example1.jpg">
   
   
  Figure 1: 3 bus example
</div>





```math
Y_{bus}=\begin{bmatrix}
-j10+(-j15) & -(-j10) & -(-j15)\\
-(-j10) & -j12+(-j10) & -(-j12)\\
-(-j15) & -(-j12) & -j15+(-j12)
\end{bmatrix}
```

```math
Y_{bus}=\begin{bmatrix}
-j25 & j10 & j15\\
j10 & -j22 & j12\\
j15 & j12 & -j27
\end{bmatrix}
```

```math
S_{2}=-2.5+j0.8
```

```math
P_{3}=2.0
```

```math
S_{i}=V_{i}\sum_{j=1}^{n} Y_{ij}^*V_{j}^* \;\;\;\;\;\; n=3
```

```math
(P_{G1}-2)+jQ_{G1}=V_{1}(j25V_{1}^*-j10V_{2}^*-j15V_{3}^*)
```

```math
-2.5+j0.8=V_{2}(j22V_{2}^*-j10V_{1}^*-j12V_{3}^*)
```

```math
2+jQ_{G3}=V_{3}(j27V_{3}^*-j12V_{2}^*-j15V_{1}^*)
```

To solve the power flow, there is no need to write equations for the bus 1(slack) but unknowns are voltage magnitude and angle at bus 2, and voltage angle at bus 3 (3 equations are enough).

```math
-2.5=\vert V_{2} \vert (1(\operatorname{Re}(j10)\cos(\theta_{2}-0)+\operatorname{Im}(j10)\sin(\theta_{2}-0)) \\

+\vert V_{2} \vert(\operatorname{Re}(-j22)\cos(\theta_{2}-\theta_{2})+\operatorname{Im}(-j22)\sin(\theta_{2}-\theta_{2})) \\

+1.1(\operatorname{Re}(j12)\cos(\theta_{2}-\theta_{3})+\operatorname{Im}(j12)\sin(\theta_{2}-\theta_{3}))) \\

=10\vert V_{2} \vert \sin(\theta_{2})+13.2\vert V_{2} \vert \sin(\theta_{2}-\theta_{3})
```

```math
2=16.5 \sin(\theta_{3})+13.2\vert V_{2} \vert \sin(\theta_{3}-\theta_{2})
```

```math
0.8=-10\vert V_{2} \vert \cos(\theta_{2})+22\vert V_{2} \vert^2 -13.2\vert V_{2} \vert \cos(\theta_{2}-\theta_{3})
```

Using the Newton-Raphson method, solution is:
```math
\theta_{2}=-4.882\degree \\
\theta_{3}=1.461\degree \\
\vert V_{2} \vert = 1.083 \\
Q_{G3}=1.967 \\
P_{G1}=2.5 \\
Q_{G1}=-2.285
```