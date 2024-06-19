# Lecture 7 - Regression, Continued

(no lecture videos)

$\newcommand{\ind}{\perp\!\!\!\!\perp}$
$\newcommand{\indicator}{1\!\!1}$

## Linear Regression and Lack of Identifiability
Consider the following model:  
$\bold{Y} = \bold{X}\beta + \varepsilon$  
with:
1. $\bold{Y} \in \R^n$ (depdendent variables), $\bold{X} \in \R^{n \times p}$ (deterministic design)
2. $\beta \in \R^p$, unkown
3. $\varepsilon \sim N_n(0, \sigma^2 I_n)$

Previously, we assumed that $\bold{X}$ had rank $p$, so we could invert $\bold{X}^\top\bold{X}$.  
What if rank $\ne p$? E.g. if $p > n$?  
Then $\beta$ is not identified, estimation of $\beta$ is vain (unless we add more structure).  

What about prediction? $\bold{X}\beta$ is still defined.  
$\hat{\bold{Y}}$: orthogonal projection of $\bold{Y}$ onto the linear span of the columns of $\bold{X}$.  
$\hat{\bold{Y}} = \bold{X}\beta = \bold{X}(\bold{X}^\top\bold{X})^{\dagger} \bold{X}^\top \bold{Y}$
- where $A^{\dagger}$ stands for the ([Moore-Penrose](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse)) pseudo-inverse of a matrix A.

Similairly as before, if $k = rank(\bold{X})$:  
- $\frac{||\hat{\bold{Y}} - \bold{Y}||^2_2}{\sigma^2} \sim \chi^2_{n-k}$
- $||\hat{\bold{Y}} - \bold{Y}||^2_2 \ind \hat{\bold{Y}}$

In particular:  
- $E[||\hat{\bold{Y}} - \bold{Y}||^2_2] = (n-k)\sigma^2$

Unbiased estimator of the variance:  
- $\hat{\sigma}^2 = \frac{1}{n-k} ||\hat{\bold{Y}} - \bold{Y}||^2_2$

## Linear Regression in High Dimension
Consider again the model:  
$\bold{Y} = \bold{X}\beta + \varepsilon$,

with:
1. $\bold{Y} \in \R^n$ (depdendent variables), $\bold{X} \in \R^{n \times p}$ (deterministic design)
2. $\beta \in \R^p$, unkown
3. $\varepsilon \sim N_n(0, \sigma^2 I_n)$

For each $i$, $X_i \in R_p$ is the vector of covariates of the $i$-th individual.  
If $p$ is too large ($p > n$), there are too many parameters to be estimated (overfitting model), although some covariates may be irrelevant.  
Solution: reduction of the dimension.  

**Idea:** Assume that only a few coordinates of $\beta$ are nonzero (but we don't know which).  
Based on the sample, select a subset of covariates and estimate the corresponding coordinates of $\beta$.  
For $S \subseteq \set{1, ..., p}$, let:  
- $\hat{\beta}_S \in \argmin_{t \in \R^s} ||\bold{Y} - \bold{X}_S \bold{t}||^2$
    - where $\bold{X}_S$ is the submatrix of $\bold{X}$ obtained by keeping only the covariates indexed in $S$.

Select a subset $S$ that minimizes the prediction error penalized by the complexity (or size) of the model:  
$||\bold{Y} - \bold{X}_S \hat{\beta}_S||^2 + \lambda ||S||$,  
where $\lambda > 0$ is the tuning paramater.  
- If $\lambda = 2\hat{\sigma}^2$, this is the *Mallow's* $C_p$ or *AIC Criterion*.  
- If $\lambda = \hat{\sigma}^2 \log n$, this is the *BIC criterion*.  

Each of these criteria is equivalent to finding $\beta \in \R^p$ that minimizes:  
$||\bold{Y} - \bold{X}\bold{b}||^2_2 + \lambda||\bold{b}||_0$,  
where $||\bold{b}||_0$ is the number of nonzero coefficients of $\bold{b}$.  

This is a computationally hard problem: nonconvex and requires to compute $2^n$ estimators (all the $\hat{\beta}_S$ for $S \subseteq \set{1, ..., p}$)

**Lasso estimator:**  
replace  
$||\bold{b}||_0 = \Sigma_{j=1}^p \indicator \set{b_j = 0}$  
with  
$||\bold{b}||_1 = \Sigma_{j=1}^p |b_j|$  
and the problem becomes convex.  
$\hat{\beta}^L \in \argmin ||\bold{Y} - \bold{X}\bold{b}||^2 + \lambda ||\bold{b}||_1$
- where $\lambda > 0$ is the tuning parameter.

How to choose $\lambda$?  
Difficult, but a good choice will lead to an estimator $\hat{\beta}$ that is very close to $\beta$ and will allow to recover the subset $S^*$ of all $j \in \set{1, ..., p}$ for which $\beta_j \ne 0$, with high probability 


## Nonparametric regression
In the linear setup, we assumed $Y_i = \bold{X}_i \beta + \varepsilon_i$, where $\bold{X}_i$ are deterministic.  

This has to be understood as working *conditionally* on the design matrix $\bold{X}$.  

And this is to assume that $E[Y_i | X_i]$ is a **linear function** of $X_i$, which is not true in general.  

Let $f(X) = E[Y_i | \bold{X}_i = x]$, $x \in \R_p$.  
How do we estimate the function $f$ ? 

**Let $p=1$ in the following:**
- One can make a parametric assumption on $f$.
    - E.g. $f(x) = a + bx$, or $f(x) = a + bx + cx^2$, or $f(x) = e^{a+bx}$, ...
- The problem reduces the estimation to a finite number of parameters.
- LSE, MLE, and all the previous theory for the linear case can be adapted.  
- But what if we do not make any such parameteric assumptions on $f$?


Assume $f$ is smooth enough: $f$ can be well approximated by a piecewise constant function.

**Idea:** local averages.  
For $x \in \R: f(t) \approx f(x)$ for $t$ close to $x$.  

For all $i$ such that $X_i$ is close enough to $x$, 
$Y_i \approx f(x) + \varepsilon_i$.  
Estimate $f(x)$ as the average of all $Y_i$'s for which $X_i$ is close enough to $x$.  

"Algorithm":
- Let $h > 0$ be the window's size (or bandwidth)  
- Let $I_x = \set{i = 1, ..., n : |X_i - x| < h}$
- Let $\hat{f}_{n,h}(x)$ be the average of $\set{Y_i: i \in I_x}$
- $\hat{f}_{n,h}(x) =$  
$\frac{1}{|I_x|}\Sigma_{i \in I_x} Y_i$ if $I_x \ne \empty$  
$0$ otherwise.
    - (the average of all y values in the set of observations with distance to $x$ less than the bandwidth)

How to choose $h$?
- If $h \to 0$: overfitting the data (too small of window, not enough generalization)
- If $h \to \infty$: underfitting (too big of window, too much generalization - if "$h = \infty$" you consider all observations which is just the sample average)
- If the smoothness of $f$ is known (i.e., quality of local approximation of $f$ by piecewise constant functions): There is a *good* choice of $h$ depending on that smoothness.  
- If smoothness of $f$ is unkown - other techniques (e.g., cross validation)  

Examples:
- k nearest neighbors - average y values of the of k closest data points to the given input x is used to estimate x.  
-  https://en.wikipedia.org/wiki/Nonparametric_regression

