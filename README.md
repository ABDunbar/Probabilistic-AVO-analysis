# Probabilistic AVO analysis
Section 4.3.16 of Quantitative Seismic Interpretation. Avseth, et al, 2005

The Jupyter notebook presents an example of probabilistic AVO analysis taken from the 2005 book, 'Quantitative Seismic Interpretation' by Avseth et.al.  The method involves calculating a reflection coefficient, R, and a gradient, G, over reflection angles 0-40 degrees.  The calculation of R and G is through an equation that takes in as input P-wave, S-wave and density for two adjacent rocks (cap-rock over reservoir).

The Vp, Vs, rho for each rock type is Monte-Carlo sampled from the probability density function derived from the histogram distribution of each log in this given well.  The random sampling quantifies the uncertainty associated with each input log and thus the AVO response for each rock pair.
![Probabilistic_AVO_analysis](https://user-images.githubusercontent.com/37248267/210221454-73b99a1b-7c90-4df7-953e-01f1ffc3f9a0.png)
