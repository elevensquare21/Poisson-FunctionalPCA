# Poisson-FunctionalPCA
Model accelerometer data under a non-homogeneous Poisson framework using multi-level function PCA technique

We leverage minute-level accelerometer data from two large randomized trials from the School of Medicine at the University of California San Diego (UCSD): MENU a randomized diet intervention among overweight women and Reach for Health a trail examining the use of metformin and a lifestyle intervention in overweight, breast cancer survivors. These two trials provide us with accelerometer data on 578 overweight women as well as a rich array of health status measurements including biomarkers and quality of life. Our model starts by assessing the underlying trends in physical activity behind the accelerometer recordings by using a non-homogeneous Poisson model. The intensity function, which parametrizes this model, is estimated using functional principal component analysis techniques tailored to our particular context. The extremely large variation of the actual data from the intensity curve motivates us to develop a formal hypothesis test for over-dispersion. As the test confirms our suspicion that a simple non-homogeneous Poisson process is insufficient, we propose two approaches to modify our original model in order to account for over-dispersion: 1. compound Poisson process 2. quasi-Poisson process.

Paper: Statistical and Computational Methods for Analyzing Accelerometer Data

File order:
1. typical day
2. start time analysis
3. minutes lost due to truncation
4. better curve registration
5. multilevel FPCA
6. plot PCs, covariance etc
7. level2 analysis
8. total variation analysis
9. menu VS rfh
10. 2 models for over-dispersion
11. L2 norm and Theta plot
12. regression analysis
13. model evaluation
14. simulation for random day
15. edgeworth simulation
