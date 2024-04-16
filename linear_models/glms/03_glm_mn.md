
# Generalized Linear Models, continued

Notes based on [the book by Nelder and McCullagh, 1989](https://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf).


## GLMs for Binary Data

M&N use "binary data" to refer to data that is based on binary (0 or 1) events, so this includes both Bernoulli and binomial data (of course, Bernoulli being a special case of the latter).

### Binary Responses
If $Y$ can only take 2 possible values, denoted 0 and 1 for convenience, we may write:  
- $pr(Y_i = 0) = 1 - \pi_i$ (probability of failure)  
- $pr(Y_i = 1) = \pi_i$ (probability of success) 

The principal objectivee of a statistical analysis is to investigate the relationship between the response probability $\pi = \pi(x)$ and the explanatory variables $x$.

Observations may be grouped in some way (e.g., "covariate classes" - where groups of observations share the same values for sets of covariates). Then, their data can be aggregated so that we record total successes and group size for each group.  
Then the responses have form $y_1/m_1, \dots, y_n/m_n$, where $m_i$ is the number of observations/trials in the $i$-th group (the "binomial index").

Note that for models applied to grouped data, asymptotic approximations for models can be based on either $m \to \infty$ or $N \to \infty$ - only the latter is appropriate for ungrouped data of course.

### Binomial Distribution

Arises when $Y$ are non-negative counts bounded above by a fixed value ($m$).

Notation $Y \sim B(m, \pi)$ - "index" $m$ and parameter $\pi$.

Suppose that $Y_1$ and $Y_2$ are independent Poisson r.v. with means $\mu_1$ and $\mu_2$. Then $Y_1 + Y_2$ is a Poisson r.v. with mean $\mu_1 + \mu_2$.  
The conditional distribution of $Y_1$ given $Y_1 + Y_2 = m$ is:  
$$
pr(Y_1 = y | Y_1 + Y_2 = m) = {m \choose y} \pi^y (1 - \pi)^{m - y}
$$
where $\pi = \mu_1 / (\mu_1 + \mu_2)$.  

The conditional distribution depends only on the Poisson means, and not on $\mu_1 + \mu_2$.

**Bernoulli distribution** is a special case where $m = 1$.  
The binomial distribution can also be thought of as the sum of $m$ independent, homogeneous Bernoulli r.v's.

The moment and cumulant generating functions of the binomial distribution are:  
- $M_y(\xi) = E \exp(\xi Y) = 1 - \pi + \pi \exp(\xi)$.  
- $K_y(\xi) = \log M_y(\xi) = \log{(1 - \pi + \pi \exp(\xi))}$

The moment and cumulant generating functions of $Y_1 + \cdots + Y_m$ are thus:  
- $(1 + \pi + \pi \exp(\xi))^m$
- $m \log{(1 - \pi + \pi \exp(\xi))}$

From the Taylor expansion of the last expression, we find the first three cumulants are:  
- $\kappa_1 = m \pi$  
    - i.e., $\mu = m \pi$
- $\kappa_2 = m \pi (1 - \pi)$
    - i.e., $V(\mu) = m \pi (1 - \pi)$
- $\kappa_3 = m \pi (1 - \pi)(1 - 2\pi)$
- *All cumulants of $Y$ have form $m$ times polynomial in $\pi$.*

Note that, in practice it is usually known only that there is variability among the $\pi$'s, and the complete set $\pi_1, \dots, \pi_m$ is unknown. Thus we often regard $\pi_1, \dots, \pi_m$ as independent r.v.'s with mean $\bar{\pi}$.  
Then, whatever the distribution of $\pi_i$, $Y_i \sim B(1, \bar{\pi})$ and $Y = Y_1 + \cdots + Y_m \sim B(m, \bar{\pi})$.


**Poisson limit**
- suppose that $\pi \to 0$ and $m \to \infty$ such that $\mu = m \pi$ remains fixed or tends to a constant.  
- Then the cumulant generating function of $Y$ tends to: 
$$
\frac{\mu}{\pi} \log{(1 + \pi(\exp(\xi) - 1))} \to \mu(\exp(\xi) - 1)
$$

This is the cumulant generating function for a Poisson r.v. with mean $\mu$.

### Logistic transformation
The *empirical logistic transformation* transforms $Y$ to attempt to achieve approximate additivity in linear logistic models.  
Suppose $Y \sim B(m, \pi)$ and that we require an approximately unbiased estimate of the log odds:  
$\lambda = \log(\frac{\pi}{1 - \pi})$  
It is natural to begin by trying transformations of the form:  
$\log \frac{Y + c}{m - Y + c}$, for some constant $c \gt 0$.  
For the particular choice $c = 1/2$, this is the *empirical logistic transformation*:  
$$
\log \frac{Y + 1/2}{m - Y + 1/2}
$$
which has property $E(Z) = \lambda + O(m^{-2})$.  

### Models for binary responses

We assume that the dependence of $\pi$ on $(x_1, \dots, x_p)$ occurs through the linear combination $\eta = X \beta$.  
We must specify a link function $g$ that maps $\pi \in (0, 1)$ to the whole real line $(-\infty, \infty)$.  

Thus we have:  
- Random component: $Y \sim B(m, \pi)$  
- Systematic component: $g(\pi_i) = \eta_i = X \beta$  
- Link function: there are four commonly used in practice:  
    1. Logit/logistic function: $$g_1(\pi) = \log(\pi / (1 - \pi))$$ (most common in practice)
        - inverse link function: $$g_1^{-1}(\eta) = \frac{\exp(\eta)}{1 + \exp(\eta)}$$ 
    2. The probit/inverse Normal function: $$g_2(\pi) = \phi^{-1}(\pi)$$, where $\phi$ is the standard Normal cdf.  
        - inverse link function: $$g_2^{-1}(\eta) = \phi(\eta)$$
    3. The complementary log-log function: $$g_3(\pi) = \log(-\log(1 - \pi))$$
        - inverse link function: $$g_3^{-1}(\eta) = 1 - \exp(-\exp(\eta))$$
    4. The log-log function: $$g_4(\pi) = -\log(-\log(\pi))$$ (seldom used)
        - inverse link function: $$g_4^{-1}(\eta) = \exp(-\exp(-\eta))$$


The logistic and probit functions are almost linearly related over interval $0.1 \le \pi \le 0.9$, so it is usually difficult to distinguish between them in a goodness-of-fit test. All asymptotic, approximate theory presented applies regardless of the choice of link.  
The logistic function is a convenient choice since it is easily interpretable as the log of the odds ratio. It also has an advantage when used with data collected retrospectively (see below).

### Parameter interpretation

If a linear logistic model is used to model binary data, e.g. with 2 covariates, we have the model:  
$$
g(\mu) = \eta \Rightarrow \log(\frac{\pi}{1 - \pi}) = \beta_0 + \beta_1 x_1
$$

Where the response is interpreted as the log odds of a success.  
The model can also be written in terms of the odds of a success:  
$$
\frac{\pi}{1 - \pi} = \exp(\beta_0 + \beta_1 x_1)
$$  

Finally, the model can be written in terms of the probability of success:
$$
\pi = \frac{\exp(\beta_0 + \beta_1 x_1)}{1 + \exp(\beta_0 + \beta_1 x_1)}
$$
- note that this is the inverse link function, $g^{-1}$.

Thus, the effect of a unit increase in $x_2$, while holding $x_1$ constant, is to:  
- increase the log odds of success by an amount $\beta_2$ (i.e., a 1 unit increase in $x_2$ increases the log odds of success by $\beta_2$)  
- increase the odds of success multiplicatively by $\exp(\beta_2)$ (i.e., a 1 unit increase in $x_2$ multiplies the odds of success by $\exp(\beta_2)$) 

On the probability scale, the effect on $\pi$ of a unit change in $x_2$ is more complicated because it depends on the values of both $x_1$ and $x_2$.  
- The derivative of $\pi$ wrt $x_2$ is $\frac{\partial \pi}{\partial x_2} = \pi(1 - \pi) \beta_2$.
- So a small change in $x_2$ has a larger effect on $\pi$ when $\pi$ is near 0.5 than when it is near 0 or 1.
- "Perhaps the simplest device to assist in the presentation of conclusions is to give the graph of $\pi(\eta) = \exp(\eta)/(1 + \exp(\eta))$ and to state the effect on $\eta$ of changes in $x_2$. The effect on the probability can then be read from the graph. This method works equally well whatever the link function used."

### Retrospective sampling

*M&N pg. 111*

Differences on the logistic scale can be estimated regardless of whether the data are samples *prospectively* or *retrospectively*.
- Prospective: data on individuals are collected as the study progresses. Retrospective: individuals are selected and then data are collected from records of past events.


### Log likelihood

The responses $y_1, \dots, y_n$ are assumed to be the observed values of independent random variables $Y_1, \dots, Y_n$, where $Y_i \sim B(m_i, \pi_i)$.  

Then the log likelihood as a function of the n-vector $\pi$ can be written as:  
$$
l(\pi;y) = \sum_{i=1}^n [y_i \log(\frac{\pi_i}{1 - \pi_i} + m_i \log(1 - \pi_i))]
$$
(ignoring $\sum \log{m_i \choose y_i}$ since it does not depend on $\pi$).

The first derivative is:  
$$
\frac{\partial l}{\partial \pi_i} = \frac{y_i - m_i \pi_i}{\pi_i (1 - \pi_i)}
$$

Using the chain rule, we can also write:  
$\frac{\partial l}{\partial \beta_r} = \sum_{i=1}^n \frac{y_i - m_i \pi_i}{\pi_i (1 - \pi_i)} \frac{\partial \pi_i}{\partial \beta_r}$
- it is convenient to express $\partial \pi_i/\partial \beta_r$ as a product:  
$\frac{\partial \pi_i}{\partial \beta_r} = \frac{d \pi_i}{d \eta_i} \frac{\partial \eta_i}{\partial \beta_r} = \frac{d \pi_i}{d \eta_i} x_{ir}$
- Thus,   
$\frac{\partial l}{\partial \beta_r} = \sum_i \frac{y_i - m_i \pi_i}{\pi_i (1 - \pi_i)} \frac{d \pi_i}{d \eta_i} x_{ir}$

The fisher information for $\beta$ is:  
$$
-E[\frac{\partial^2 l}{\partial \beta_r \partial \beta_s}] = \sum_i \frac{m_i}{\pi_i (1 - \pi_i)} \frac{\partial \pi_i}{\partial \\beta_r} \frac{\partial \pi_i}{\partial \beta_s} \\
= \sum_i m_i \frac{(d \pi_i / d \eta_i)^2}{\pi_i(1 - \pi_i)} x_{ir} x_{is} \\
= (X^T W X)_{rs}
$$

where $W$ is a diagnoal weight matrix given by:
$$
W = diag (m_i (\frac{d \pi_i}{d \eta_i})^2 / \pi_i (1 - \pi_i))
$$
- note that for the logistic link, $\partial l / \partial \beta = X^\top (Y - \mu)$ (when written in matrix notation) and $W = diag(m_i \pi_i (1 - \pi_i))$.  

#### Parameter estimates

Using the Newton-Raphson method, the parameter estimates are obtained by:  
- Given initial estimates $\hat{\beta}_0$, compute $\hat{\pi}_0$ and $\hat{\eta}_0$.  
- Then define the adjusted dependent variable $Z$ with components:  
$z_i = \hat{\eta}_i + \frac{y_i - m_i \hat{\pi}_i}{m_i} \frac{d \eta_i}{d \pi_i}$  
- ML estimates satisfy:  
$X^\top W X \hat{\beta} = X^\top W Z$
- This can be solved iteratively using least squares methods, so the revised estimate is:  
$\hat{\beta}_1 = (X^\top W X)^{-1} X^\top W Z$

## Deviance function

The residual deviance is defines as twice the difference between the maximum achievable log likelihood and the log likelihood of the fitted model:  
$$
D(y; \hat{\pi}) = 2l(\tilde{\pi}; y) - 2l(\hat{\pi}; y) \\ 
= 2 \sum_i [y_u \log(y_i / \hat{\mu_i}) + (m_i - y_i) \log(\frac{m_i - y_i}{m_i - \hat{\mu_i}})]
$$



## Bias and Precision of Estimates

To a first order approximation, ML estimates are unbiased with asymptotic variance equal to the inverse Fisher information matrix.  

For large $n$, $E(\hat{\beta} - \beta) = O(n^{-1})$ and $cov(\hat{\beta}) = (X^\top W X)^{-1} \{1 + O(n^{-1})\}$.  
These approx results are also true for the limit where $n$ is fixed and $m \to \infty$.  

## Sparseness

When a sizeable proportion of observed counts $m_i$ are small, we say the data is sparse.  
Sparseness can affect the deviance function and Pearson's statistic, which will fail to have the properties required of a goodness-of-fit statistic.  
(Note that if at least $n$ is large, the asymptotic approximations for the ML estimates are still valid.)  


For example, if $Y_i \sim B(1, \pi_i)$ (bernoulli), and a logistic model has been fit yielding  
$\hat{\pi}_i = \exp(x_i^\top \beta) / (1 + \exp(x_i^\top \beta))$,  
then the residual deviance function is:  
$D = 2 \sum_i [y_i \log(y_i / \hat{\pi}_i) + (1 - y_i) \log((1 - y_i) / (1 - \hat{\pi}_i)) - y_i \log(\hat{\pi}_i / (1 - \hat{\pi}_i)) - \log(1 - \hat{\pi}_i)$.  
But since $y_i$ is either 0 or 1, $y \log y = (1 - y) \log(1 - y) = 0$, and $\log(\hat{\pi}_i / (1 - \hat{\pi}_i)) = x_i^\top \hat{\beta}$, the deviance function simplifies to:  
$$
D = -2 \hat{\beta}^\top X^\top Y - 2 \sum_i \log(1 - \hat{\pi_i})
$$

Thus, $D$ is a function of $\hat{\beta}$, and conditional on $\hat{\beta}$, $D$ has a degenerate distribution (i.e., it is essentially deterministic - only exactly degenerate for $m = 1$).


## Over-dispersion

Over-dispersion refers to when the variance of the response $Y$ exceeds the nominal binomial variance, $m \pi (1 - \pi)$.  
Over-dispersion is likely to be common in practice. Depending on the application, the binomial variance may just be a small component of the total variance in $Y$.  

Over-dispersion can arise in many ways, but a common cause is clustering in the population (e.g., families, schools, neighborhoods, etc.).  

Assume for simplicity that cluster size $k$ is fixed and that the $m$ sampled individuals come from $m/k$ clusters.  
For the $i$-th cluster, the number of positive responses is $Z_i \sim B(k, \pi_i)$.  
Then the total number of positive responses is  
$Y_i = Z_1 + Z_2 + \dots + Z_{m/k}$.  

Writing $E(\pi_i) = \pi$ and $var(\pi_i) = \tau^2 \pi (1 - \pi)$, it may be shown that the unconditional mean and variance of $Y$ are:  
- $E(Y) = m\pi$
- $var(Y) = m \pi (1 - \pi) (1 + (k - 1)\tau^2)$
    - Rewriting, $var(Y) = \sigma^2 m \pi (1 - \pi)$, where dispersion parameter $\sigma^2 = 1 + (k - 1)\tau^2$.
    - Note that the dispersion parameter depends on cluster size and on the variablity of $\pi$ from cluster to cluster.  

This derivation forces $\sigma^2$ to be $1 \le \sigma^2 \le k \le m$ since $0 \le \tau^2 \le 1$.  

To accomodate under-dispersion, it is often desirable to extend this so that $\sigma^2$ can be between 0 and 1 as well.

The beta-binomial distribution is sometimes used as an alternative to model over-dispersion, but it is generally unwise to rely on a specific form of over-dispersion where the assumed form is chosen mostly for mathematical convenience (*M&N pg 126*).

Over-dispersion can only occur if $m \gt 1$, otherwise the mean determines all cumulants (dispersion parameter is 1 for a bernoulli model).  


### Parameter estimation

Let the effect of over-dispersion be to inflate the variance by $\sigma^2$ (i.e., $var(Y) = \sigma^2 m \pi (1 - \pi)$).

Then such a model can be fit using the same methods as for the binomial model, and we will continue to assume that the binomial distribution applies.  

The main difference:  
- The covariance matrix of $\hat{\beta}$ is replaced with: $$cov(\hat{\beta}) \approxeq \sigma^2 (X^\top W X)^{-1}$$
    - where $W = diag(m_i (\frac{d \pi_i}{d \eta_i})^2 / \pi_i (1 - \pi_i))$, from the paramater estimation section (4.4.2), which also gets calculated during IRLS iterations.  


We must also **estimate the dispersion paramater** $\sigma^2$ (or sometimes referred to generally as $\phi$), which is analogous to estimating the $\sigma^2$ parameter in a Normal-theory linear model. 

In the absence of replication (discussed below), or if the number of degrees of freedom for replication is small, an estimate of $\sigma^2$ may be based on a weighted version of the residual sum of squares.  
If the fitted model is correct, the following is approximately unbiased for $\sigma^2$ provided that $p$ is small compared to $n$:  
$$
\tilde{\sigma}^2 = \frac{1}{n - p} \frac{(y_i - m_i \hat{\pi}_i)^2}{m_i \hat{\pi}_i (1 - \hat{\pi}_i)} = X^2 / (n - p)
$$
(i.e., Pearson's statistic divided by the residual degrees of freedom).  

The estimated covariance matrix of $\hat{\beta}$ is then:  
$$
\widehat{var}(\hat{\beta}) = \tilde{\sigma}^2 (X^\top W X)^{-1}
$$

#### Replication
If there is replication, then for each covariate value $x$, several observations $(y_1, m_1), \dots, (y_r, m_r)$ are made, which are independent and essentially identically distributed, except that the indicies $m_i$ may be unequal.  

Then the estimate of $\pi$ based on this single covariate class is $\tilde{\pi} = y_. / m_.$.  
The expected value of the within-class weighted sum of squares is:  
$\sum_{j=1}^r (y_j - m_j \tilde{\pi})^2 / m_j = (r - 1)\sigma^2 \pi (1 - \pi)$,  
so that:
$$
s^2 = \frac{1}{r - 1} \sum_j \frac{(y_j - m_j \tilde{\pi})^2}{m_j \tilde{\pi}(1 - \tilde{\pi})}
$$

These covariate class estimators are then pooled together to obtain the replication estimate of dispersion on $\sum(r - 1)$ degrees of freedom.











