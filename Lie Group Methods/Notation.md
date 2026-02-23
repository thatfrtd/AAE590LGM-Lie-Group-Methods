2D Column Vector ```latex \Vtwo{a}{b}```: $\newcommand{\Vtwo}[2]{\begin{pmatrix}#1 \\ #2 \end{pmatrix}} \Vtwo{a}{b}$

$2\times2$ Matrix ```latex \Mtwo{a}{b}{c}{d}```: $\newcommand{\Mtwo}[4]{\begin{pmatrix}#1 & #2 \\ #3 & #4 \end{pmatrix}}\Mtwo{a}{b}{c}{d}$

Symmetric $2\times2$ Matrix ```latex \symMtwo{a}{b}{c}```: $\newcommand{\symMtwo}[3]{\Mtwo{#1}{#3}{#3}{#2}}\symMtwo{a}{b}{c}$

Diagonal $2\times2$ Matrix ```latex \diagMtwo[1]{x}{y}```: $\newcommand{\diagMtwo}[3][0]{\symMtwo{#2}{#3}{#1}}\diagMtwo[1]{x}{y}$
	Default off-diagonal is zero: $\diagMtwo{x}{y}$

Constant Diagonal $2\times2$ Matrix ```latex \cdiagMtwo[1]{x}```: $\newcommand{\cdiagMtwo}[2][0]{\diagMtwo[#1]{#2}{#2}}\cdiagMtwo[1]{x}$
	Default off-diagonal is zero: $\cdiagMtwo{x}$

3D Column Vector ```latex \Vthree{a}{b}{c}```: $\newcommand{\Vthree}[3]{\begin{pmatrix}#1 \\ #2 \\ #3 \end{pmatrix}} \Vthree{a}{b}{c}$

$3\times3$ Matrix ```latex \Mthree{a}{b}{c}{d}{e}{f}{g}{h}{i}```: $\newcommand{\Mthree}[9]{\begin{pmatrix}#1 & #2 & #3 \\ #4 & #5 & #6 \\ #7 & #8 & #9 \end{pmatrix}}\Mthree{a}{b}{c}{d}{e}{f}{g}{h}{i}$

Symmetric $3\times3$ Matrix ```latex \symMthree{a}{b}{c}{d}{e}{f}```: $\newcommand{\symMthree}[6]{\Mthree{#1}{#4}{#6}{#4}{#2}{#5}{#6}{#5}{#3}}\symMthree{a}{b}{c}{d}{e}{f}$

Diagonal $3\times3$ Matrix ```latex \diagMthree[1]{x}{y}{z}```: $\newcommand{\diagMthree}[4][0]{\symMthree{#2}{#3}{#4}{#1}{#1}{#1}}$
$\diagMthree[1]{x}{y}{z}$
	Default off-diagonal is zero: $\diagMthree{x}{y}{z}$

Constant Diagonal $3\times3$ Matrix ```latex \cdiagMthree[1]{x}```: $\newcommand{\cdiagMthree}[2][0]{\diagMthree[#1]{#2}{#2}{#2}}\cdiagMthree[1]{x}$
	Default off-diagonal is zero: $\cdiagMthree{x}$

Adjoint ```latex \Ad{X}{x}```: $\newcommand{\Ad}[2]{\text{Ad}_{#1}(#2)}\Ad{X}{x}$
adjoint ```latex \ad{\xi}```: $\newcommand{\ad}[1]{\text{ad}(#1)}\ad{\xi}$
Exponential map ```latex \Exp{\xi}```: $\newcommand{\Exp}[1]{\text{Exp}(#1)}\Exp{\xi}$
Logarithmic map ```latex \Log{\X}```: $\newcommand{\Log}[1]{\text{Log}(#1)}\Log{X}$

Differential ```latex \d```: $\newcommand{\d}{\text{d}}\d$ 

Derivative ```latex \der[2]{x}{y}```: $\newcommand{\der}[3][]{\frac{\d^{#1} #2}{\d #3^{#1}}}\der[2]{x}{y}$
	Default power is one: $\der{x}{y}$