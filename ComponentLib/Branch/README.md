# Branch Model

**Note: Branch model not yet implemented**

Transmission lines and different types of transformers (traditional, Load Tap-Changing transformers (LTC) and Phase Angle Regulators (PARs)) can be modeled with a common branch model. 
The most common circuit that is used to represent the transmission line model is $`\pi`$ circuit as shown in Figure 1 along with the flow directions (Sending end = FROM bus and Receiving end = TO bus).



<div align="center">
   <img align="center" src="Documentation/Figures/TL.jpg">
   
   
  Figure 1: Transmission line $`\pi`$ equivalent circuit
</div>






``` math
Z'=R+jX
```
``` math
Y'=G+jB
```


where $`R`$ is line series resistance, $`X`$ is line series reactance, $`B`$ is line shunt charging, and $`G`$ is line shunt conductance. As can be seen from Figure 1 total $`B`$ and $`G`$ are separated between two buses. 
Following equation can be written:

```math
\begin{bmatrix}
I_{S}\\
-I_{R}
\end{bmatrix}=Y_{TL}\begin{bmatrix}
V_{S}\\
V_{R}
\end{bmatrix}
```

where:

```math
Y_{TL}=\begin{bmatrix}
\dfrac{1}{R+jX}+\dfrac{G+jB}{2} & -\dfrac{1}{R+jX}\\
-\dfrac{1}{R+jX} & \dfrac{1}{R+jX}+\dfrac{G+jB}{2}
\end{bmatrix}
```

The branch model can be created by adding the ideal transformer in series with the $`\pi`$ circuit as shown in Figure 2 where $`\tau`$ is a tap ratio magnitude and $`\theta_{shift}`$is the phase shift angle.

<div align="center">
   <img align="center" src="Documentation/Figures/branch.jpg">
   
   
  Figure 2: Branch equivalent circuit
</div>


The branch admitance matrix is then:

```math
Y_{BR}=\begin{bmatrix}
(\dfrac{1}{R+jX}+\dfrac{G+jB}{2})\dfrac{1}{\tau^2} & -\dfrac{1}{R+jX}\dfrac{1}{\tau e^{-j\theta_{shift}}}\\
-\dfrac{1}{R+jX}\dfrac{1}{\tau e^{j\theta_{shift}}} & \dfrac{1}{R+jX}+\dfrac{G+jB}{2}
\end{bmatrix}
```
## Branch Contribution to Buses i and j

Branch b between buses $`i`$ and $`j`$ contributes to 4 elements of the system $`Y`$ matrix, $`Y_{ii}`$, $`Y_{ij}`$, $`Y_{ji}`$, and $`Y_{jj}`$.
Using the following formulations:

```math
P_{i}= \vert V_{i} \vert \sum_{j=1}^{n}\vert V_{j} \vert (G_{ij}\cos\theta_{ij}+B_{ij}\sin\theta_{ij})
```

```math
Q_{i}= \vert V_{i} \vert \sum_{j=1}^{n} \vert V_{j} \vert (G_{ij}\sin\theta_{ij}-B_{ij}\cos\theta_{ij})
```

branch contributions are as follows:
```math
P_{i}= \vert V_{i} \vert (\vert V_{j} \vert (G_{ij}\cos\theta_{ij}+B_{ij}\sin\theta_{ij})+\vert V_{i} \vert (G_{ii}\cos\theta_{ii}+B_{ii}\sin\theta_{ii}))
```
```math
Q_{i}= \vert V_{i} \vert (\vert V_{j} \vert (G_{ij}\sin\theta_{ij}-B_{ij}\cos\theta_{ij})+\vert V_{i} \vert (G_{ii}\sin\theta_{ii}-B_{ii}\cos\theta_{ii})
```
```math
P_{j}= \vert V_{j} \vert (\vert V_{i} \vert (G_{ji}\cos\theta_{ji}+B_{ji}\sin\theta_{ji})+\vert V_{j} \vert (G_{jj}\cos\theta_{jj}+B_{jj}\sin\theta_{jj}))
```
```math
Q_{j}= \vert V_{j} \vert (\vert V_{i} \vert (G_{ji}\sin\theta_{ji}-B_{ji}\cos\theta_{ji})+\vert V_{j} \vert (G_{ji}\sin\theta_{jj}-B_{jj}\cos\theta_{jj}))
```

where $`B_{ii}`$, $`B_{ij}`$, $`B_{ji}`$,$`B_{jj}`$, and $`G_{ii}`$, $`G_{ij}`$, $`G_{ji}`$,$`G_{jj}`$ can be calculated from the $`Y_{BR}`$ elements as follows:

``` math
G_{ii}+jB_{ii}=Y^{BR}_{ii}
```
``` math
G_{ij}+jB_{ij}=Y^{BR}_{ij}
```
``` math
G_{ji}+jB_{ji}=Y^{BR}_{ji}
```
``` math
G_{jj}+jB_{jj}=Y^{BR}_{jj}
```
