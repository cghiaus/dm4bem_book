# Dynamic Models for Building Energy Management

by Christian Ghiaus

- Researcher ID: ORCID [0000-0001-5561-1245](https://orcid.org/0000-0001-5561-1245)
- SciProfiles [2970335](https://sciprofiles.com/profile/2970335)
- Scopus [6603390490](https://www.scopus.com/authid/detail.uri?authorId=6603390490)
- Web of Science [K-1307-2012](https://www.webofscience.com/wos/author/record/1651371)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem_book/blob/main/LICENSE)

State-space representation is widely used for developing control algorithms. This Jupyter Book shows how thermal transfer in buildings can be modeled by complex thermal networks, assembled from elementary thermal networks, that are converted in state-space representation ([Ghiaus 2013](https://hal.archives-ouvertes.fr/hal-03605823/document), [Ghiaus 2021](https://doi.org/10.1007/978-3-030-76477-7_5)). These steps are implemented by using `dm4bem` module written in _Python 3.9_ (and tested on _Pyhton 3.11_). The book does not cover the development of control algorithms.


__Quick overview__

The workflow is presented in the Jupyter Notebook on [inputs and simulation](tutorials/pd05simulation.ipynb).


__Prerequisites__

It is assumed that readers have a foundational knowledge at the undergraduate level in the areas of linear algebra ([Strang, G. 2023](https://math.mit.edu/~gs/linearalgebra/ila6/indexila6.html)), heat transfer ([RE2020 2021](https://rt-re-batiment.developpement-durable.gouv.fr/IMG/pdf/annexeiv_arrete_4_aout_2021.pdf)), and Python programming ([Docs.Pyhton 2024](https://docs.python.org/3/tutorial/index.html)).


__Notations used for values of quantities__

This book uses the writing conventions for SI unit symbols and names recommanded by the *International Bureau of Weights and Measures* ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true), [Gőbel et al. 2006](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3)) and *National Institute of Standards and Technology* ([Thomson and Taylor. 2008](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf)).


__Reproducibility__

The Jupyter Notebooks can be run interactively on [mybinder.org](https://mybinder.org) by pushing the __launch binder__ button:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem_book/HEAD)

__Contents__

```{tableofcontents}
```

__References__

1. [Ghiaus, C. (2013)](https://doi.org/10.1016/j.energy.2012.10.024). Causality issue in the heat balance method for calculating the design heating and cooling loads, *Energy* 50: 292-301, [hal-03605823](https://hal.archives-ouvertes.fr/hal-03605823/document)

2. [Ghiaus, C. (2021)](https://doi.org/10.1007/978-3-030-76477-7_5). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198, [HAL 03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)

3. [BIPM (2019)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true). The International System of Units (SI), 9th edition, licence CC-BY-3.0

4. [Gőbel, E., Mills, I., Wallard,  A. (2006)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3). A concise summary of the International System of Units, the SI

5. [Thomson, A., Taylor, B. N. (2008)](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf). Guide for the use of the international System of units (NIST Special Publication 811․ 2008 Edition). National Institute of Standards and Technology, US Government Printing Office.

6. [Strang, G. (2023)](https://math.mit.edu/~gs/linearalgebra/ila6/indexila6.html). Introduction to Linear Algebra, 6th ed., ISBN 978-17331466-7-8

7. [RE2020 (2021)](https://rt-re-batiment.developpement-durable.gouv.fr/IMG/pdf/annexeiv_arrete_4_aout_2021.pdf). Annexe IV : Règles « Th-Bat 2020 » - données d’entrée au calcul de la performance énergétique

8. [Docs.Pyhton (2024)](https://docs.python.org/3/tutorial/index.html) The Python Tutorial