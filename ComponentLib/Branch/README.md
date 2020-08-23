# Branch Model

Transmission lines and different types of transformers (traditional, Load Tap-Changing transformers (LTC) and Phase Angle Regulators (PARs)) can be modeled with a common branch model. 
The most common circuit that is used to represent the transmission line model is $`\pi`$ circuit as shown in Figure 1 along with the flow directions (Sending end = FROM bus and Receiving end = TO bus).



<div align="center">
   <img align="center" src="https://gitlab.pnnl.gov/gridkit/gridkit/-/raw/my_develop_new/Documentation/Figures/TL.jpg">
   
   
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
   <img align="center" src="https://gitlab.pnnl.gov/gridkit/gridkit/-/raw/my_develop_new/Documentation/Figures/branch.jpg">
   
   
  Figure 2: Branch equivalent circuit
</div>


The branch admitance matrix is then:

```math
Y_{BR}=\begin{bmatrix}
(\dfrac{1}{R+jX}+\dfrac{G+jB}{2})\dfrac{1}{\tau^2} & -\dfrac{1}{R+jX}\dfrac{1}{\tau e^{-j\theta_{shift}}}\\
-\dfrac{1}{R+jX}\dfrac{1}{\tau e^{j\theta_{shift}}} & \dfrac{1}{R+jX}+\dfrac{G+jB}{2}
\end{bmatrix}
```