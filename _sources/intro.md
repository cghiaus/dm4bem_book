# Dynamic Models for Building Energy Management

by Christian Ghiaus
- ORCID [0000-0001-5561-1245](https://orcid.org/0000-0001-5561-1245)
- SciProfiles [2970335](https://sciprofiles.com/profile/2970335)
- Scopus [6603390490](https://www.scopus.com/authid/detail.uri?authorId=6603390490)
- Web of Science [K-1307-2012](https://www.webofscience.com/wos/author/record/1651371)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem_book/blob/main/LICENSE)

Developing control algorithms requires dynamic models of the processes. This Jupyter Book explores the modeling of thermal transfer in buildings through thermal networks. It outlines the construction of elementary thermal networks and their assembly to create the thermal network of complex systems, such as a building. Additionally, it presents the transformation of these thermal networks into state-space representation, which can then serve for simulations, as demonstrated here.

This book does not cover the development of control algorithms.

The book uses `dm4bem` module written in _Python 3.9_ (and tested on _Pyhton 3.11_).


__Quick overview__

The workflow is presented in the Jupyter Notebook on [inputs and simulation](tutorials/pd05simulation.ipynb).


__Prerequisites__

It is assumed that readers have a foundational knowledge at the undergraduate level in the areas of linear algebra, heat transfer, and Python programming.


__Notations used for values of quantities__

This book uses the writing conventions for SI unit symbols and names recommanded by the *International Bureau of Weights and Measures* ([BIPM 2019](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true), [Gőbel et al. 2006](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3)) and *National Institute of Standards and Technology* ([Thomson A. et al. 2008](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf)).


__Reproducibility__

The Jupyter Notebooks can be run interactively on [mybinder.org](mybinder.org) by pushing the __launch binder__ button:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem_book/HEAD)

__Contents__

```{tableofcontents}
```

__References__
1. [BIPM (2019)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true) The International System of Units (SI), 9th edition, licence CC-BY-3.0

2. [Gőbel, E., Mills, I., Wallard,  A. (2006)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3). A concise summary of the International System of Units, the SI

3. [Thomson, A., Taylor, B. N. (2008)](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf). Guide for the use of the international System of units (NIST Special Publication 811․ 2008 Edition). National Institute of Standards and Technology, US Government Printing Office.