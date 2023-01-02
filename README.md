# Probabilistic AVO analysis
Section 4.3.16 of Quantitative Seismic Interpretation. Avseth, et al, 2005

The Jupyter notebook presents an example of probabilistic AVO analysis taken from the 2005 book, 'Quantitative Seismic Interpretation' by Avseth et.al.  The method involves calculating a reflection coefficient, R, and a gradient, G, over reflection angles 0-40 degrees.  The calculation of R and G is through an AVO approximation equation that takes in as input P-wave, S-wave and density for two adjacent rocks (cap-rock over reservoir).

The Vp, Vs, rho correlated triplet for each rock type is randomly sampled from the probability density function derived from the histogram distribution of each log in this given well.  The repeated random sampling (according to the pdf) quantifies the uncertainty associated with each input log and thus the AVO response for each rock pair.

![Probabilistic_AVO_analysis](https://user-images.githubusercontent.com/37248267/210221454-73b99a1b-7c90-4df7-953e-01f1ffc3f9a0.png)

A note about the code.  The code is written in Python and is my own except where I give reference (e.g. 2015 Alessandro Amato del Monte, for the well plot).  Avseth et.al. do provide some Matlab code but without a Matlab license or experience, I haven't used it. Therefore, there may be glaring mathematical and/or coding blunders that I hope you won't blame me for because you used it without your own quality control.  A link to the text book (which is fantastic and highly recommended) is here...
https://www.cambridge.org/core/books/quantitative-seismic-interpretation/EB6A36B78CCF07187723F6F5364EDCF8
