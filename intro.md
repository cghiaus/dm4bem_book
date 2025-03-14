# Dynamic Models for Building Energy Management

by Christian Ghiaus (Researcher ID: [ORCID](https://orcid.org/0000-0001-5561-1245), [SciProfiles](https://sciprofiles.com/profile/2970335), [Scopus](https://www.scopus.com/authid/detail.uri?authorId=6603390490), [Web of Science](https://www.webofscience.com/wos/author/record/1651371), [HAL](https://cv.hal.science/cghiaus))

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem_book/HEAD)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem_book/blob/main/LICENSE)
![Hits](https://hits.sh/https://cghiaus.github.io/dm4bem_book.svg)


> “Il y a trois sortes de savoir : le [savoir](https://fr.m.wikipedia.org/wiki/Savoir) proprement dit, le [savoir-faire](https://fr.m.wikipedia.org/wiki/Savoir-faire) et le [savoir-vivre](https://fr.m.wikipedia.org/wiki/Civilité) ; les deux derniers dispensent assez bien du premier." (_There are three types of knowledge: [knowledge](https://en.m.wikipedia.org/wiki/Knowledge) itself, [know-how](https://en.m.wikipedia.org/wiki/Procedural_knowledge), and [know how to live](https://en.m.wikipedia.org/wiki/Etiquette) (or [soft skills](https://en.m.wikipedia.org/wiki/Soft_skills)); the last two quite adequately dispense with the first_), [Talleyrand (1754-1838)](https://en.m.wikipedia.org/wiki/Charles_Maurice_de_Talleyrand-Périgord)


State-space representation is widely used for developing control algorithms. This Jupyter Book shows how thermal transfer in buildings can be modeled by assembling thermal networks that are then converted to state-space representation  ([Ghiaus 2013](https://hal.archives-ouvertes.fr/hal-03605823/document), [Ghiaus 2021](https://hal.science/hal-03578578/document)). These steps are implemented by using [dm4bem](tutorials/dm4bem.py) module written in _Python 3.8_ and tested on _Python 3.11_. The book does not cover the development of control algorithms.


__Quick overview__

The specifics of this book are:

- Modeling heat transfers (conduction, convection, short wave and long wave radiation, advection) through a matrix representation of thermal circuits.
- Formulating the thermal load calculation as a control problem.
- Obtaining thermal circuits of walls through spatial discretization.
- Assembling thermal circuits.
- Transforming thermal circuits into state-space representations.
- Providing examples for nonlinear models and control algorithms.

The workflow is presented in [Assembling thermal circuits](tutorials/pdREADME.md) and an example is given in the Jupyter Notebook on  [inputs and simulation](tutorials/pd05simulation.ipynb) (see section _Python script_). [A short example on GitHub](https://github.com/cghiaus/dm4bem_toy_model) shows the workflow.


__Prerequisites__

It is assumed that readers have a foundational knowledge at the undergraduate level in linear algebra ([Strang, G. 2023](https://math.mit.edu/~gs/linearalgebra/ila6/indexila6.html)), heat transfer ([Incropera et al. 2007](https://hyominsite.files.wordpress.com/2015/03/fundamentals-of-heat-and-mass-transfer-6th-edition.pdf), [RE2020 2021](https://rt-re-batiment.developpement-durable.gouv.fr/IMG/pdf/annexeiv_arrete_4_aout_2021.pdf)), and Python programming ([Docs.Python 2024](https://docs.python.org/3/tutorial/index.html)).


__Notations used for values of quantities__

This book uses the writing conventions for SI unit symbols and names recommanded by the *International Bureau of Weights and Measures* ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true), [Gőbel et al. 2006](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3)) and *National Institute of Standards and Technology* ([Thomson and Taylor. 2008](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf)).

Some rules:
- [_Unit symbols_](https://en.m.wikipedia.org/wiki/Unit_of_measurement) are in roman type and [_quantity symbols_](https://en.m.wikipedia.org/wiki/Physical_quantity) are in italic; the unit symbol is placed after the numerical value and a space is left between the numerical value and the unit symbol ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true) §2.1, pp.129-142), e.g. $h = 10 \, \mathrm{W \, m^{−2} \, K^{−1}}$ or $h = 10\, \mathrm{W·m^{−2}·K^{−1}}$ or $h$ = 10 W/(m²·K).
- Symbols for units formed from other units by multiplication are indicated by means of either a half-high (that is, centered) dot or a space, e.g. W/(m⋅K) or W/(m K).
- A prefix symbol attached to a unit symbol constitutes a new inseparable symbol, forming a multiple or submultiple of the unit concerned ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true), §3, pp.143-144), e.g. $1 \, \mathrm{mK} = 10^{-3} \, \mathrm{K}$ while $1 \, \mathrm{m·K} = 1 \, \mathrm{m \ K} = 1 \, \mathrm{m} · 1 \, \mathrm{K}.$
- When writing the value of a quantity as the product of a numerical value and a unit, both the number and the unit may be treated by the ordinary rules of algebra ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true), §5.4, pp.148-151, [Thomson and Taylor, 2008](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf), §7.1 p.15), e.g.:
    - In $T = 273.15 \, \mathrm{K}$, the number $273.15 = T /\mathrm{K}$ is the numerical value of thermodynamic temperature $T.$
    - The numerical value of a temperature expressed in degrees Celsius, $\theta$, is related to the numerical value of the thermodynamic temperature expressed in kelvins, $T$, by the relation $\theta /\mathrm{°C} = T/ \mathrm{K} − 273.15$ or $\theta$/(°C) = $T$/(K) − 273.15.
    - The ordinate of a graph or the heading of a table is labeled $Temperature,$ $T/(10³$ K $)$, where $T$ is the thermodynamic temperature and K is the unit symbol for kelvin. If the ordinate value of a point on a curve of the graph or the entry in a table is 0.273, then the corresponding temperature is $T/(10³$ K $)$ = 0.273 or $T$ = 273 K.
    - The logarithm ordinate of a graph or the heading of a table is labeled $\log_{10}(T /\textrm{(K)})$. If the ordinate value of a point on a curve of the graph or the entry in a table is 2.436..., then $\log_{10}(T/\textrm{(K)}) = 2.436...$ (or $T/\textrm{(K)} = 10^{2.436...}$ ) and the corresponding temperature is $T$ = 273 K.
    - Symbol % (percent) is used for the number 0.01, e.g., the emmisivity is $\varepsilon$ = 0.85 = 85 %.

__Nomenclature__

| Symbol  |                 |
|:--------|:----------------|
| ±       | Plus-minus sign |
|⁻ ¹ ² ³ ⁴| Superscripts    |
| · ×     | Multiplication  |

|Symbol| SI unit |          | Quantity |
|:-----|:--------|----------|:----------|
| _t_  | s       |          | Time      |
| _w_  | m       |          | Widh      |
| _θ_  | °C      |          | Temperature |
| _q_  | W       |          | Heat flow rate |
|_A, S_| m²      |          | Surface area |
| _V_  | m³      |          | Volume |
| _T_  | K       |          | Thermodynamic temperature|
| _T_  | °C      |          | Temperature source|
|_Q̇, Φ_| W       |          | Heat flow rate source |
|      |         |          |                |
| _λ_  | W/(m·K) |W·m⁻¹·K⁻¹ | Thermal conductivity|
| _ρ_  | kg/m³   |kg·m⁻³    | Density |
| _c_  | J/(kg·K)|J·kg⁻¹·K⁻¹| Specific heat capacity |
| _h_  | W/(m²·K)|W·m⁻²·K⁻¹ | Heat transfer coefficient |
| _C_  | J/K     |J·K⁻¹     | Thermal capacity |
| _E_  | W/m²    |W·m⁻²     | Irradiance |
| _G_  | W/K     |W·K⁻¹     | Thermal conductance |
| _R_  | K/W     |K·W⁻¹     | Thermal resistance |
| _ṁ_  | kg/s    |kg·s⁻¹    | Mass flow rate |
| _V̇_  | m³/s    |m³·s⁻¹    | Volumetric flow rate |

__Reproducibility__

The results presented are [reproducible](https://en.m.wikipedia.org/wiki/Reproducibility). The Jupyter Notebooks can be run interactively on [mybinder.org](https://mybinder.org) by pushing the __launch binder__ button:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem_book/HEAD)

__Licences__

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem_book/blob/main/LICENSE)

The text of this book is under Creative Commons Attribution 4.0 International (CC BY 4.0).

The software is under MIT Licence.


__GitHub repository:__ https://github.com/cghiaus/dm4bem_book


__Contents__

```{tableofcontents}
```

__References__

1. [Ghiaus, C. (2013)](https://doi.org/10.1016/j.energy.2012.10.024). Causality issue in the heat balance method for calculating the design heating and cooling loads, *Energy* 50: 292-301, [hal-03605823](https://hal.archives-ouvertes.fr/hal-03605823/document)

2. [Ghiaus, C. (2021)](https://doi.org/10.1007/978-3-030-76477-7_5). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198, [hal-03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)

3. [BIPM (2019)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true). The International System of Units (SI), 9th edition, licence CC-BY-3.0

4. [Gőbel, E., Mills, I., Wallard,  A. (2006)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3). A concise summary of the International System of Units, the SI

5. [Thomson, A., Taylor, B. N. (2008)](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf). Guide for the use of the international System of units (NIST Special Publication 811․ 2008 Edition). National Institute of Standards and Technology, US Government Printing Office.

6. [Strang, G. (2023)](https://math.mit.edu/~gs/linearalgebra/ila6/indexila6.html). Introduction to Linear Algebra, 6th ed., ISBN 978-17331466-7-8

7. [Incropera, F. P., DeWitt, D. P., Bergman, T. L., Lavine, A. S. (2007)](https://hyominsite.files.wordpress.com/2015/03/fundamentals-of-heat-and-mass-transfer-6th-edition.pdf). Fundamentals of Heat and Mass Transfer, 6th Edition. John Wiley.

8. [RE2020 (2021)](https://rt-re-batiment.developpement-durable.gouv.fr/IMG/pdf/annexeiv_arrete_4_aout_2021.pdf). Annexe IV : Règles « Th-Bat 2020 » - données d’entrée au calcul de la performance énergétique

9. [Docs.Python (2024)](https://docs.python.org/3/tutorial/index.html) The Python Tutorial

