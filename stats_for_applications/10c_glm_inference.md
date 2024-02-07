# More on GLMs

Following along to [Heather Turner's GLM course](https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf) as a road map though  [the book by Nelder and McCullagh, 1989](https://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf).



Starts as mostly a review of GLM theory notes in 10a, but expanded upon a bit to bring in more general treatment of the dispersion parameter and prior weights.

## Exponential Family


Recall the general exponential family's PDF follows the form:  
$$
f_\theta(y) = h(y) \exp(\eta(\theta)T(y) - B(\theta))
$$
where $h(y)$ and $B(\theta)$ are known functions (see 10a notes for explanation - $h$ is a change of measure and $B$ is normalizing).

For example,  
Consider $X \sim N(\mu, \sigma^2)$, s.t. $\theta = (\mu, \sigma^2)$, which is a 2-parameter exponential family (unless $\sigma$ is known, in which case it is a 1-parameter).  
- The density is:  
$f_\theta(y) = \exp(\frac{\mu}{\sigma^2}y - \frac{1}{2\sigma^2}y^2 - \frac{\mu^2}{2 \sigma^2}) \frac{1}{\sigma \sqrt{2 \pi}}$  
Which forms the 2-paramater exponential family with:  
    - $\eta_1 = \frac{\mu}{\sigma^2}$, $\eta_2 = \frac{-1}{2\sigma^2}$  
    - $T_1(y) = y$, $T_2(y) = y^2$,  
    - $B(\theta) = \frac{\mu^2}{2\sigma^2} + \log(\sigma \sqrt{2 \pi})$  
    - $h(y) = 1$  


The exponential family can be generalized to the overdispersed exponential family, which has a structure that is nice for statistical inference.  
$$
f_\theta(y) = h^*(\phi, y) \exp(\frac{\eta(\theta)T(y) - B(\theta)}{a(\phi)})
$$
For some specific functions $a, b, c$.  

If $\eta$ is identity than the distiribution is said to be in canonical or natural form. Further, if $T$ is identity and $\phi$ is known this is said to be a "one-parameter exponential family". If $\phi$ is not know, this may or may not be a two-parameter exponential family.
$$
f_\theta(y) = \exp(~~ \frac{y \theta - b(\theta)}{a(\phi)} + c(y, \phi) ~~)
$$


As shown in previous notes (10a), we can derive $E(Y)$ (or, more commonly, $E(Y|X)$) and the $Var(Y)$ (or $Var(Y|X)$) from $b(\theta)$, the "cumulant generating function" (a result related to Bertlett's identities).
$$
E(Y) = b'(\theta) = \mu
$$
$$
Var(Y) = a(\phi) b''(\theta) = a(\phi) V(\mu)
$$
- Where $V(\mu)$ is known as the **variance function**. $b''(\theta)$ depends only on the canonical parameter, and hence the mean, as expressed by $V(\mu)$.

- Typically, $a(\phi) = \frac{\phi}{w}$, where $w$ is the provided weight on the response variable, typically just $1$, though included here for generality, and $\phi$ is the dispersion parameter (sometimes written $\sigma^2$). $a(\phi)$ is independent of $\theta$, depending only on $\phi$.

    - The variance function is especially important because it determines which specific family of distributions withing the exponential family a distribution belongs. For example, for the Normal distribution, $V(\mu) = 1$ s.t. $Var(Y) = \phi$ (usually written $\sigma^2$). For the Poisson distribution $V(\mu) = \mu$ and for the Binomial distribution $V(\mu) = \mu(1 - \mu)$, and $\phi$ is typically taken to be $1$ for these distributions which gives the Poisson and Binomial $Var(Y)$'s that we all know.
    - Thus, the variance is a function of the mean in a specific way. But it can only depend on the covariates through the mean $\mu$ that they imply. 

## GLM
Think of a GLM as modelling $\mu = E[Y|X]$ by $g(\mu) = X \beta$. The GLM is often characterized by 3 components:  
1. A random/stochastic component: this refers to the distribution from which $Y|X$ is drawn - a distribution with expected value $\mu$ and potentially a scale/dispersion parameter as well. We assume independence and constant variance typically.
2. A systematic component: this refers to to $\eta = X \beta$, the linear predictor.  
3. The link $g$ between the random and systematic components: $g(\mu) = \eta$


## Canonical Links

First note that in the world of GLMs, we often write $X_i^\top \beta = \eta$, introducing a new symbol to describe the linear predictor.  

Recall that the link function is $g$ in the GLM $g(\mu_i) = X_i^\top \beta$, where $X_i, \beta \in \R^p$, s.t. $g$ maps the mean/location parameter $\mu$ from its support to $\R$, so that it can be modelled by $X_i^\top \beta \in R$.  

Recall that the **canonical link** is defined as:  
$$
g(\mu_i) = (b')^{-1}(\mu_i) = \theta_i = X_i^\top \beta = \eta
$$

The canonical link provides us with useful statistical properties but we do not need to use it. If we do not use the canonical link, then we consider $\theta$ to be a function of $X$, recall:  
$\theta_i(X_i) = (b')^{-1}(\mu_i(X_i)) = (b')^{-1}(g^{-1}(X_i^\top \beta)) := h(X_i^\top \beta) = h(\eta)$  
(thus, when $g = (b')^{-1}$ as in the canonical case, $h$ is identity).



## Estimation of Parameters

Recall we estimate $\theta$ by maximizing the log likelihood function:  
$$
\ell_n(\beta; \bold{Y}, \bold{X}) = \sum_i \log f_\beta(Y, X) = \sum_{i=1}^n \frac{Y_i \theta_i - b_i(\theta)}{a(\phi)} \\
$$
$$
= \sum_{i=1}^n \frac{Y_i h(X_i^\top \beta) - b(h(X_i^\top \beta))}{a(\phi)}
$$
- (note we drop $c(y_i, \phi_i)$ because it does not depend on $\beta$, and thus doesn't matter for optimization)
- note that sometimes this is written for a single observation and so the summation term is dropped


The MLE for $\beta_j$ is obtained by solving $(\nabla \ell_n)_j  = 0$ (the score equations):  
- $\frac{\partial \ell_n}{\partial \beta_j} = \sum_{i=1}^n \frac{Y_i - \mu_i}{a(\phi_i) V(\mu_i)} \frac{X_{ij}}{g'(\mu_i)} = 0$
- note that $1 / g'(\mu_i) = 1 / (d \mu / d \eta) = d \eta / d \mu$, i.e., $\mu'(\eta)$  

- Also note that "$\phi_i$ disappears" (Nelder and McCullagh pg 42) - if $a(\phi_i) = \phi_i / w_i$, then $\phi_i$ is constant across observations and can be dropped but $w_i$ needs to be added into the above equation since it varies across observations (?) unless $w$ is 1 in which case it can essentially be ignored. Leaving:  
$$\frac{\partial \ell_n}{\partial \beta_j} = \sum_{i=1}^n \frac{w_i (Y_i - \mu_i) }{V(\mu_i)} X_{ij} \mu'(\eta) = 0$$



In the previous notes (10b) we saw that the Fisher-scoring method can be used given the gradient and Hessian of the log likelihood, and that this can itself be rewritten as a weighted least squares problem given the current estimates at iteration $k$, $\beta^{(k)}$, where $Z$ is the response:   
$$
Z_i^{(k)} = \eta + \frac{(Y_i - \mu_i^{(k)})}{\mu'(\eta)}
$$

And $W$ is the weight matrix:  
- $W^{(k)} =  diag \frac{\mu'(\eta)^2}{Var(\mu) a(\phi_i)}$
- However, as discussed above $a(\phi_i) = \phi_i / w_i$, and $\phi_i$ drops off so we are left with:  
$$W^{(k)} = diag \frac{w_i ~ \mu'(\eta)^2}{V(\mu)}$$

Recall, the algorithm is:  
1. Fix $\beta^{(k)}$ and set $\mu_i^{(k)} = g^{-1}(X_i^\top \beta^{(k)})$
2. Calculate the adjusted dependent responses: $Z_i^{(k)}$.  
3. Compute the weights $\bold{W}^{(k)} = \bold{W}(\beta^{(k)})$.  
4. Regress (i.e., WLS) $\bold{Z}^{(k)}$ on the design matrix $\bold{X}$ with weight $\bold{W}^{(k)}$ to derive a new estimate:  
$$\beta^{(k + 1)} = (\bold{X}^\top \bold{W}^{(k)} \bold{X})^{-1} \bold{X}^\top \bold{W^{(k)}} Z^{(k)}$$  
5. Repeat until convergence.  

To perform this algorithm we just need:
- Data ($X$ and $Y$), initial values for $\beta$, and any weights $w$.
- $\mu'(\eta)$, i.e., $d(g^{-1})/d \eta$ - so the derivative of the link function.
- $V(\mu)$ - the variance as a function of the mean for the exponential family






## Measuring Goodness of Fit

### Deviance 

**Motivation**  

(Nelder and McCullagh Chapter 2.1.2)  
When we estimate the parameters of a model, we want to define a measure of goodness of fit between the observations and the model's fitted values so that we can assess the precision of the estimates, and hence the validity of our model.  
With GLMs we are primarily interested in the **log likelihood** function $l(\mu; y)$ (in these notes also sometimes written $\ell_n(\beta; Y, X)$ in the regression context).  
$$
l(\mu; y) = \sum_i \log f_i(y_i; \theta_i)
$$
- the sum of the individual logged probability density functions, where $\mu$ and $y$ are vectors (and if only 1 obs. the summation term can be dropped)

As seen already, we fit GLMs by maximizing the log likelihood...

But instead of using the log likelihood as a goodness-of-fit criterion, there are advantages to using a particular linear function of log likelihoods called the **scaled deviance**:  
$$
D^*(y, \mu) = 2 l(y; y) - 2l(\mu; y)
$$

- For exponential family models, $l(\mu = y; y)$ is the maximum likelihood achievable for an exact fit, where the fitted values are equal to the observed ($\hat{\mu} = y$).  
- $l(\mu, y)$ is the log likelihood of the chosen model, as we are used to thinking about it.  


Since $l(y,y)$ does not depend on the parameters of the model, if we were to minimize $D^*(y, \mu)$ this term could be dropped (it's effectively a constant) which means minimizing $D^*$ would be equivalent to maximizing the log likelihood $l(\mu, y)$.  

As an example, consider linear regression / Gaussian GLM with known variance:  
- We have the probability density
$$
f(y; \mu) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp(- \frac{(y - \mu)^2}{2 \sigma^2})
$$  
- The log likelihood
$$
l(\mu; y) = -\frac{1}{2} \log(2 \pi \sigma^2) - (y - \mu)^2 /(2 \sigma^2)
$$
- But setting $\mu = y$ yields the following, which is the maximum achievable log likelihood:  
$$
l(y; y) = -\frac{1}{2} \log(2 \pi \sigma^2)
$$
- Plugging these into $D^*$ we get the scaled deviance is:
$$
D^*(y; \mu) = 2(l(y;y) - l(\mu;y)) = (y - \mu)^2 / \sigma^2
$$
- Up to the constant $\sigma^2$, this is the Residual Sum of Squares and this implies that minimizing $D^*$ is equivalent to minimizing the RSS, i.e. least squares.



**Theory**  
(Nelder and McCullagh Chaper 2.3)


Let $\hat{\mu}$ be the vector of fitted values from the GLM, so that $Y - \hat{\mu}$ is the residual.  
Measures of discrepancy or goodness of fit between $\hat{\mu}$ and $y$ can be formed in various ways, but with GLMs we tend to consider the *deviance* - the logarithm of a ratio of likelihoods (for the reasons shown above).  

The **null model** is a GLM with one parameter (i.e., just an intercept), which represents a common $\hat{\mu}$ for all the observations. Since no covariates are incorporated into the model, this model attributes all the variation in the $Y$'s to the random component.  

The **full model** or saturated model has $n$ parameters (one per observation), and thus the $\hat{\mu}$'s derived from it match the data ($Y$) exactly, so that all of the variation in $Y$ is attributed to the systematic component.

In GLM theory, we define the discrepancy of a fit, i.e., the **deviance**, to be proportional to 2 times the difference between the maximum log likelihood achievable and that achieved by the model under investigation (i.e., the log likelihood of the full model - where the fitted values for $\hat{\mu}$ are $Y$ - minus the log likelihood of the fitted model).  
$$
D^*(Y, \hat{\mu}) = D(Y, \hat{\mu}) / \phi = 2 (\ell_{full} - \ell_{model})
$$
$$
= \frac{2w}{\phi} \sum_{i=1}^n (Y_i [\theta(Y_i) - \theta(\hat{\mu}_i)]- b(\theta(Y_i)) + b(\theta(\hat{\mu}_i))) 
$$
- given $a(\phi_i) = \phi_i / w_i$, and note that $D^*$ is the "scaled deviance" and $D$ is just the "deviance".
- as usual, $b$ is the cumulant generating function for the given exponential family
- $\theta$ is the canonical parameter. If the canonical link is used, than $\theta = X \beta$ and it doesn't need to be thought of as a function of the data. Otherwise, $\theta = (g \circ b')^{-1}(X \beta)$ 


Some examples of $D$:  
- Normal: $w \sum (Y - \hat{\mu})^2$, the sum of squared residuals
- Poisson: $2w \sum (y \log(y / \hat{\mu}) - (y - \hat{\mu}))$
- Binomial: $2w \sum(y \log(y / \hat{\mu}) + (m - y) \log[(m - y)/m - \hat{\mu}])$
- Gamma: $-2w \sum( \log(y / \hat{\mu}) + (y - \hat{\mu}) / \hat{\mu})$
- Inverse Gamma: $w\sum (y - \hat{\mu})^2 / (\hat{\mu}^2 y)$


Note: it is a common misconception that $D^* \sim \chi^2_{n - p}$. This is only true in some special circumstances, but in general the $\chi^2$ approximations of the deviance are not good even as $n \to \infty$ (Nelder and McCullagh, pg 36).
The misconception may be due to the likelihood ratio test, see that section.  

However, we can utilize $\chi^2$ approximations of the *differences* between deviances for nested models.

**Analysis of Deviance**  

This is related to the idea of the likelihood ratio test, but not exactly the same.  

The idea is to compare the null model (or any nested model) with the proposed model. If the proposed model's deviance is much smaller than the null model, it suggests that adding terms to the model improved it.  
This is what Nelder and McCullagh mean when they say we can look at the difference between the deviances of nested models. Doing so with the null model deviance and that of the proposed model, the difference is approximately distributed as $\chi^2_{p-1}$ ($p-1$ d.o.f. because we have d.o.f. null - d.o.f. proposed = (n - 1) - (n - p) = p - 1, assuming p includes the intercept). 

But this can be done with any two nested models:  
Say $M_2$ has $p_2$ predictors and $M_1$ has $p_1 < p_2$ predictors, so that $M_1$ is **nested** within $M_2$ (for example, $M_1$ may be the null model and $M_2$ a model with some covariates, or $M_1$ may have some covariates and $M_2$ has those but more).  
We can test the null hypothesis that the *extra* coefficients of $M_2$ are simultaneously zero. If we can't reject this null, this suggests the additional covariates offer no benefit when included.  
$H_0: \beta_{p_1 + 1} = \dots = \beta_{p_2} = 0$  
$H_1: \beta_j \ne 0$, for any $\beta_j \in \beta_{p_1 + 1}, \dots, \beta_{p_2}$

The test statistic is:  
$T_n = D^*_{p_1} - D^*_{p_2} \sim \chi^2_{p_2 - p_1}$  

So, for the $\alpha$ % test, we reject the null if $T_n$ is greater than the $1 - \alpha$ quantile of the $\chi^2_{p_2 - p_1}$ distribution.  

*However*, as seen above $D^*$ depends on $\phi$, which we have to estimate, and thus we cannot actually rely on $D^*_{p_1} - D^*_{p_2}$ being $\chi^2_{p_2 - p_2}$.  

Instead, we compute the $F$ test statistic:  
$F = \frac{(D^*_{p_1} - D^*_{p_2}) / (p_2 - p_1)}{D^*_{p_2} / (n - p_2 - 1)} = \frac{(D_{p_1} - D_{p_2})/(p_2 - p_1)}{D_{p_2} / (n - p_2 - 1)} \sim F_{p_2 - p_1, n - p_2 - 1}$

- So the dispersions $\phi$ in the $D^*$ cancel out.

(Notes pull from N&M and https://bookdown.org/egarpor/PM-UC3M/glm-deviance.html)



### Residuals

We can define multiple forms of generalized residuals:  

The **Pearson residual**:  
$$
r_P = \frac{Y - \hat{\mu}}{\sqrt{V(\mu)}}
$$
- i.e., the raw residual scaled by the estimated std. deviation


The **deviance residual**:  
$$
r_D = sign (Y - \hat{\mu}) \sqrt{d_i}
$$
- Where $D = \sum_{i=1}^n d_i$, i.e., the (unscaled) deviance is broken into a sum of values that come from each data point, so $d_i$ is the deviance to due observation $i$.  
- $sign$ indicates that the direction of the residual comes from the difference between the actual and fitted values (positive if you underestimate, negative if you overestimate)
- But the magnitue is defined just using the individual deviance. Recall the deviance formula from above is a sum across the $n$ observations, so we could just pluck out the $i$-th observation in the deviance vector before summing.
-  Note that this implies that $D = \sum_{i=1}^n (r_D)^2$ - the sume of squared deviance residuals (recall in the deviance section we already showed that the deviance for a Normal model is the SSR!)


Sometimes you will see the "response" residual, which is just $Y - \hat{\mu}$


## Estimating the (Nuisance) Dispersion Parameter

In the models discussed here we have assumed $\phi$ is a constant, although it may be unkown, and hence is a [nuisance parameter](https://en.wikipedia.org/wiki/Nuisance_parameter), because it is not directly a quantity of interest but is needed for hypothesis testing of $\beta$.  


Recall that $Var(Y) = a(\phi) b''(\theta) = \frac{\phi}{w} V(\mu)$ so $\phi = \frac{w Var(Y)}{V(\mu)}$.  


If $\phi$ is unknown, an estimate is required.
There are practical difficulties in estimating the dispersion $\phi$ by maximum likelihood and it is usually estimated by method of moments. If $\beta$ was known an unbiased estimate would be:   
$\phi = \frac{1}{n} \sum_{i=1}^n \frac{w_i (y_i - \mu_i)^2}{V(\mu_i)}$

Since $\beta$ must be estimated we obtain:  
$$
\hat{\phi} = \frac{1}{n - p} \sum_{i=1}^n \frac{w_i (y_i - \mu_i)^2}{V(\mu_i)}
$$
- where $w$ is the vector of prior weights, $\mu = g^{-1}(X\beta)$ corresponds to the final fitted values, and $n - p$ is the residual degrees of freedom.

- Note that in the Gaussian case (with weights = 1) this simplifies to:  
$\hat{\phi} = \frac{1}{n-p} \sum_{i=1}^n (y_i - \mu_i)^2$, which is equivalent to the $\hat{\sigma}^2$ formula solved for in linear regression.

Recall the generalized Pearson residuals from above:  
$r_P = \frac{Y - \hat{\mu}}{\sqrt{V(\mu)}}$, which means we could calculate dispersion via:  
$$
\hat{\phi} = \frac{1}{n - p} \sum_{i=1}^n w_i  (r_P)^2 
$$
which is how R does it in `summary.glm`.  

However, R calculates the pearson residuals as:  
where the residuals $r_P = (Y - \mu) / \mu'(\eta)$ 



## Covariance Matrix and Standard Errors for $\beta$

The GLM estimates $\hat{\beta}$ have the usual properties of maximum likelihood estimators, particularly they are asymptotically normal with covariance matrix equal to the inverse of the Fisher information matrix:  
$$\hat{\beta} \sim N(\beta, \phi(X^\top W X)^{-1}) \equiv \hat{\beta} \sim N(\beta, I^{-1})$$  

Where $I(\beta) = \phi^{-1} (X^\top W X)$ and we can use $\hat{W} = W^{(k_{last})}$, the final weights from the IRLS procedure, to estimate $W$.
- Recall that (for a single obs.) $W = diag \frac{w ~ \mu'(\eta)^2}{V(\mu)}$, and that $I = - E[H_{\ell_n}]$
- Earlier we saw that the gradient of the log likelihood (for a single obs.) was such that $\frac{\partial \ell_n}{\partial \beta_i} = \frac{(Y - \mu)}{V(\mu) a(\phi)} \mu'(\eta) X_i$  
- The second derivative takes the form: $\frac{\partial^2 \ell_n}{\partial \beta_j \partial \beta_k} = \frac{\partial^2 \ell_n}{\partial \eta^2} X_i X_j$  
- And thus the negative expectation of the second derivative is:  
$-E[H(\ell_n)_{ij}] = \frac{\mu'(\eta)^2}{V(\mu) a(\phi)} X_i X_j = \frac{w ~\mu'(\eta)^2}{V(\mu) \phi} X_i X_j$  (since $a(\phi) = \phi / w$), but we can rewrite this in terms of $W$:  
$= \frac{1}{\phi} W X_i X_j$
- Or in matrix form: $-E[H_{\ell_n}] = \frac{1}{\phi} X^\top W X$  
- The inverse, i.e. the Fisher information, is $\phi (X^\top W X)^{-1}$
- (from [Nelder and Wedderburn, 1972](https://www.medicine.mcgill.ca/epidemiology/hanley/bios601/Likelihood/NelderWedderburn1972.pdf), pg 373)

So, $$\hat{cov}(\hat{\beta}) = \phi(X^\top W X)^{-1}$$

Thus, **standard errors** for for the $\beta_j$ coefficient may be calculated as the square root of the diagonal elements of $\hat{cov}(\hat{\beta})$.  

If $\phi$ is unkown, an estimate is required.




## Wald Tests & T-Tests

We can again use the fact that (if $\phi$ is known)   
$$\hat{\beta} \sim N(\beta, \phi(X^\top W X)^{-1})$$  

And use a z-test to test the significance of a coefficient.  
Specifically, we test  
$H_0 : \beta_j = 0$ vs.   
$H_1 : \beta_j \ne 0$

with test stat:  
$z_j = \frac{\hat{\beta}_j}{\sqrt{\phi (X^\top W X)^{-1}_{jj}}}$

which is asymptotically $N(0, 1)$ under the null.  
So our test is  $\psi = \set{ |z_j| \ge q_{\alpha/2}}$ where $q_{\alpha/2}$ is the $1 - \alpha / 2$ quantile the $N(0,1)$, for the $\alpha$ % level test.   
We can calculate a p-value as $2 * (1 - \Phi(|T_n|))$, the probability of observing a more extreme test stat if the null is true.  

However, for exponential families without a know dispersion (i.e., not binomial or Poisson for whom dispersion is assumed to be $1$, but when we must use a $\hat{\phi}$), we can do a t-test instead of Wald's.  
In which case, the test statistic is the same but it is asymptotically distributed as a Student's T with $n - p$ degrees of freedom.  
(As an example, consider the T test for (Normal) linear regression)


## Confidence and Prediction Intervals

### Linear Regression 

https://online.stat.psu.edu/stat501/lesson/3/3.3

https://stackoverflow.com/questions/38109501/how-does-predict-lm-compute-confidence-interval-and-prediction-interval/38110406#38110406


**Confidence Interval**  

The confidence interval for observation $h$ generally takes the form: $\hat{y}_h \pm t_{1 - \alpha, n - 2} SE(\hat{y}_h)$, where $t_{1 - \alpha/2, n - p}$ is the $1 - \alpha/2$ quantile of the Student's T distribution with $n - p$ degrees of freedom, which corresponds to the $\alpha$ % level test. Thus, we need to calculate $SE(\hat{y}_h)$. 
- note, the T distribution has $n-p$ d.o.f. because our estimate for the variance of y, $\hat{\sigma}^2$, has the $n-p$ in the denominator (where $\hat{\sigma}^2$ is sometimes called mean squared error or MSE in this context) 

Recall that $cov(\hat{\beta}) = \sigma^2 (X^\top W X)^{-1}$, the inverse of the Fisher information matrix.  

For a single observation $h$, $\hat{y}_h = x_h^\top \hat{\beta}$.    
Thus, $var(\hat{y}_h) = var(x_h^\top \hat{\beta}) = x_h^\top cov(\hat{\beta})x_h = x_h^\top (\sigma^2 (X^\top W X)^{-1})x_h = \sigma^2 x_h^\top (X^\top W X)^{-1} x_h$  
We also estimate $\sigma^2$ such that $\hat{\sigma}^2 = \sum_i (y_i - x_i^\top \hat{\beta})^2 / (n - p)$, leaving us with the following formula:  

$SE(\hat{y}) = diag \sqrt{\hat{\sigma}^2 x_h^\top (X^\top W X)^{-1}x_h}$

For a single observation, this can be shown to be:  
$SE(\hat{y}_h) = \sqrt{\hat{\sigma}^2 (\frac{1}{n} + \frac{(x_h - \bar{x})^2}{\sum_i (x_i - \bar{x})^2})}$


**Prediction Interval**  
Now we are interested in the standard error of ${y}_{new}$, a prediction made from $x_{h(new)}$, a new observation.  
Given a new observation, we would like to be able to predict a value for the response variable ($\hat{y}_{h(new)}$) and give some estimate of uncertainty. For example, if the mean and standard deviation of $y$ at $x_{h(new)}$ were known, we could say that about 95% of the time we run this experiment the true value $y_{new}$ should like within the interval $(\mu - 2 \sigma, \mu + 2 \sigma)$.  
However, we do not actually know $\mu$ at $x_{h(new)}$ and so instead we estimate it using our fitted model, which will give us a $\hat{y}_{h(new)}$. However, we pay a cost by estimating $\mu$ with $\hat{y}_{h(new)}$ - the variance of $\hat{y}_{h(new)}$. Since different samples will lead to different parameters in the model and hence to different $\hat{y}_{h(new)}$'s, we have to take into account the variance of $\hat{y}_{h(new)}$'s sampling distribution.  
Thus, **the variance of the prediction has 2 components**:  
- The variation of $y$, i.e. $\sigma^2$, which we estimate as $\hat{\sigma}^2$ (since $y \sim N(X\beta, \sigma^2 I)$) 
- The variation due to estimating the mean $\mu$ with $\hat{y}_h$ - i.e., the variance of the sampling distribution of $\hat{y}$, which is the same SE that we see in the confidence interval.  

So for a single observation we have:  
$SE(\hat{y}_{h(new)}) = \sqrt{\hat{\sigma}^2(1 + \frac{1}{n} + \frac{(x_{h(new)} - \bar{x})^2}{\sum_i (x_i - \bar{x})^2})}$

Comparing to the confidence interval formula above, we can see that the formulas are very similair. The prediction interval has an extra $\hat{\sigma}^2$ term in the square root, since we are accounting for the variance in $y_{new}$ itself in addition to its sampling variance.

i.e.,  
$SE = X \Sigma X^\top$, where for the confidence interval $\Sigma = \Sigma_{\hat{\beta}} = \sigma^2 (X^\top X)^{-1}$, and for the prediction interval, $\Sigma = \Sigma_{\hat{\beta}} + \sigma^2 I$ (adding in the prediction uncertainty)

### GLMs

https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

It is evident that we cannot simply add and substract 2 times the SE of the reponse from the fitted value in the context of GLMs, because it may result in outcomes that violate the parametric assumptions imposed by picking a specific family of probability distributions. For example, imagine a binomial GLM with a fitted value close to 0. Subtracting 2 times the SE may result in a lower-bound that is negative, which is not allowed for the the parameter $p$ of the binomial distribution.  

A simple method is to compute the confidence bounds on the "link" scale, instead of the "response" scale (i.e., $g(\mu) \pm 2 SE(g(\hat{\mu}))$).  
Since the link function maps $\mu$ to the real line, the confidence interval on the link scale is allowed to take any value. Then we apply the inverse of the link function to map it back to the response scale.  

As we saw above with classic linear regression, since $g(\hat{\mu}) = X \hat{\beta}$,  
$$Var(g(\hat{\mu})) = Var(X \hat{\beta}) = X cov(\hat{\beta}) X^\top$$
where $cov(\hat{\beta}) = \phi(X^\top W X)^{-1}$ (inverse of Fisher information).

Then we can do: $g(\hat{\mu}) \pm  t_{1 - \alpha/2, n-p} ~ SE(g(\hat{\mu}))$ to calculate the lower and upper bounds, and then map the bounds back to response scale via $g^{-1}$

- Note, can it be shown that $Var(\hat{\mu}) = g^{-1}(Var(X \beta))$? Or what is the proof that this is a valid confidence interval?  
- I think this will help: https://stats.stackexchange.com/questions/503814/what-is-the-best-way-to-form-a-confidence-interval-for-a-weighted-sum-of-expecte


This is more or less what R does, calculating the variance of each fitted value, on the link scale, as $diag \sqrt{X S X^\top}$, where $S = \phi (X^\top W X)^{-1}$, the covariance matrix of $\hat{\beta}$.  

However, if you specify type = "response", [it multiplies this SE of the fitted value](https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/predict.glm.R#L60) by $\mu'(\eta)$... why?

- The idea is [apparently](https://stackoverflow.com/questions/44534864/what-standard-errors-are-returned-with-predict-glm-type-response-se-fi) to give a first order approximation of the true standard error via the delta method. **But since the inverse link function $g^{-1}$ is generally non-linear, the confidence interval should not be symmetric around fitted values, so getting the confidence interval on the link scale and mapping it back to the response scale (as discussed above) is the recommended way to go.**

- To convert the standard error of the fitted mean from the link to the response scale, R uses a calculation based on a [delta method](https://en.wikipedia.org/wiki/Delta_method) approximation.  
- This is what (I believe) is happening:  
    - Consider $Var(\hat{\mu}_i) = Var(g^{-1}(X_i\hat{\beta})) = Var(g^{-1}(\hat{\eta}_i)) = Var(\mu(\hat{\eta}_i))$.  
    - This is where the delta method comes in. We can apprximate $\mu(\hat{\eta})$ using first order approximation (taylor series expansion):  
    $\mu(\hat{\eta}_i) \approx \mu(\hat{\eta}_i) + \mu'(\eta)(\hat{\eta}_i - \eta)$  
    - Essentially, [this means](https://www.stata.com/support/faqs/statistics/delta-method/) that $Var(\mu(\hat{\eta}_i)) \approx \mu'(\eta)^2 Var(\hat{\eta}_i)$, because ...  
    - using the [mean value theorem](http://users.uoa.gr/~npapadat/TaylorTheoremAmStat2013.pdf):  
    $\mu(\hat{\eta}_i) - \mu(\eta) \approx \mu'(\eta)(\hat{\eta}_i - \eta)$.  
    - Note that the Central Limit Theorem tells us:  
    $\sqrt{n}(\mu(\hat{\eta}_i) - \mu(\eta)) \xrightarrow[]{D} N(0, \mu'(\eta)^\top \Sigma_{\hat{\eta}} \mu'(\eta))$, where $\Sigma_{\hat{\eta}}$ is the covariance matrix of the sampling distribution of $\hat{\eta}$ (i.e., the square of the SE of $\hat{\eta}$).  
    - Thus, the variance of the sampling distribution of $\hat{\mu}_i$ is $Var(\hat{\mu}_i) \approx \mu'(\eta)^2 \Sigma_{\hat{\eta}}$, (since we compute the standard error of the fit for a single observation, it simplifies).  
    - The standard error of the fit (i.e., the standard deviation of the sampling distribution) is the square root: $SE(\hat{\mu}_i) \approx \mu'(\eta) \sqrt{\Sigma_\eta}$.  
    The square root of the covariance matrix of the $\hat{\eta}$ we already know - it's the sampling variance calculated above, since $g(\hat{\mu}) = \hat{\eta}$.  
    This could explain why R calculates the standard error of the fit for the response scale as:  
    $SE(g(\hat{\mu})) \times \mu'(\hat{\eta})$
- Note that using this adjusted SE of the fit for confidence intervals still doesn't make sense, because (like mentioned above), the link function, and hence inverse link function, is often non-linear which means confidence intervals will be non-linear. Hence, just calculate them on the link scale and then map back to response scale.



**Prediction intervals:**  
Not as straight-forward with GLMs it seems...

One approach is bootstrapping, of course.  

It seems that another popular approach is to use quantiles from the "posterior distribution" - in other words, transform the fitted $\mu$ and $\phi$ (if applicable) to scale of the response distribution (i.e., via the inverse link) and use them to obtain the quantiles corresponding to some upper and lower probability (say, .025 and .975 for 95% coverage).  
E.g., if for $X = 2$ your estimate/fitted value for the Poisson paramater is $\hat{\lambda} = 3.4$, then you could do `qpois(p = c(0.025, 0.975), lambda = 3.4)`

See:
- https://fromthebottomoftheheap.net/2017/05/01/glm-prediction-intervals-i/

- https://cran.r-project.org/web/packages/trending/vignettes/prediction_intervals.html

- And how SAS calculates them: ["Quantile-based prediction intervals for generalized linear models"](https://support.sas.com/kb/69/692.html)


## Akiake Information Criterion

The AIC is a measure of fit that penalizes for the number of parameters $p$.  
A smaller value indicates better fit.  
$AIC = -2 \ell_{model} + 2p$  
where $\ell_{model}$ is the log likelihood of the fitted model.  
AIC can be used to compare 2 models (not necessarily nested models, for that we may use the likelihood ratio).

Since a lower AIC is better, we can see that larger log-likelihood values lead to a lower AIC (so it rewards goodness of fit, since the likelihood function measures the "likelihood" of the parameter(s) given the data), while penalizing over-fitting via adding too many predictors (hence the $+2p$).  
Thus, we can use AIC as a relative statistic to compare models.





part 2: binary data
https://statmath.wu.ac.at/courses/heather_turner/glmCourse_002.pdf
- also: confidence intervals,
- likelihood ratio test
- model interpretation


part 3: count data https://statmath.wu.ac.at/courses/heather_turner/glmCourse_003.pdf
- also: offsets,
- overdispersion & quasi-likelihood




[GLM Model Diagnostics](https://bookdown.org/egarpor/PM-UC3M/glm-diagnostics.html)

[GLM's with shrinkage](https://bookdown.org/egarpor/PM-UC3M/glm-shrink.html)  

[What important ideas came since Nelder and McCullagh's book Generalized Linear Models (a 40 year old book)?](https://stats.stackexchange.com/questions/573586/what-important-ideas-came-since-nelder-and-mccullaghs-book-generalized-linear-m)
