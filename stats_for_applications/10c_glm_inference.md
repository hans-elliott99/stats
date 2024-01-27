## Extra

### Standard Errors
From [Heather Turner GLM course](https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf):  

**Covariance Matrix and SEs**  
The GLM estimates $\hat{\beta}$ have the usual properties of maximum likelihood estimators, particularly they are asymptotically normal with covariance matrix equal to the inverse of the Fisher information matrix:  
$\hat{\beta} \sim N(\beta, \phi(X^\top W X)^{-1})$  
$\equiv \hat{\beta} \sim N(\beta, I^{-1})$  
Where $I(\beta) = \phi^{-1} (X^\top W X)$ and we can use $\hat{W} = W^{(k_{last})}$, the final weights from the IRLS procedure, to estimate $W$.

So, $\hat{cov}(\hat{\beta}) = \phi(X^\top W X)^{-1}$

Thus, standard errors for for the $\beta_j$ coefficient may be calculated as the square root of the diagonal elements of $\hat{cov}(\hat{\beta})$.  

**Dispersion**  
Note that $Var(Y) = \phi b''(\theta) = \phi V(\mu)$, so $\phi = \frac{Var(Y)}{V(\mu)}$, where we see $V$ is the variance function that describes how $Y$'s variance depends on the mean.  
(For a Gaussian model, $\phi = \sigma^2$ and $V(\mu) = 1$)

We assume that $\phi_i$ (dispersion for subject $i$) is $\phi_i = \phi / w_i$, where $w_i$ are known *prior weights* (often 1).  
Therefore, $\phi_i = \frac{w_i Var(Y)}{V(\mu_i)}$


If $\phi$ is unknown, an estimate is required.
There are practical difficulties in estimating the dispersion $\phi$ by maximum likelihood and it is usually estimated by method of moments. If $\beta$ was known an unbiased estimate would be:   
$\phi_i = \frac{1}{n} \sum_{i=1}^n \frac{w_i (y_i - \mu_i)^2}{V(\mu_i)}$

Since $\beta$ must be estimated we obtain:  
$$
\hat{\phi_i} = \frac{1}{n - p} \sum_{i=1}^n \frac{w_i (y_i - \mu_i)^2}{V(\mu_i)}
$$
- where $w$ is the vector of prior weights, $\mu = g^{-1}(X\beta)$ corresponds to the final fitted values, and $n -p$ is the residual degrees of freedom.

Note that in the Gaussian case (with weights = 1) this simplifies to:  
$\hat{\phi} = \frac{1}{n-p} \sum_{i=1}^n (y_i - \mu_i)^2$, which is equivalent to the $\hat{\sigma}^2$ formula solved for in linear regression.

But this is how R does it for all families:  
$\hat{\phi} = \frac{\sum_{i=1}^n w_i (y_i - \mu_i)^2}{n - p}$

(ie, it just ignores the $V(\mu)$ term from the above derivation...)


Need to read more.  
- These slides: https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf
- https://stats.stackexchange.com/questions/224302/how-does-r-function-summary-glm-calculate-the-covariance-matrix-for-glm-model