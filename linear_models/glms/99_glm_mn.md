
# Generalized Linear Models, continued

Notes based on [the book by Nelder and McCullagh, 1989](https://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf).


Some random-ish bits from the book that are interesting and apply broadly across different types of statistical models.


## Types of Covariates
### Continuous Covariates

*M&N pg. 51*  

"Models containing only continuous covariates are often called regression models, to be contrasted with analysis-of-variance models, which have only terms involving qualitative factors."

On functions of continuous covariates:
- In the context of GLMs, linearity refers to linearity of $\eta$ in terms of its parameters.  
- A continuous covariate $x$ in a model may be replaced by an arbitrary function $g(x)$ without destroying the linearity of the model. For example, we may add $x^2$ and $x^3$
    - *When reading about GAMs: after the transformation of a covariate into basis functions, we are essentially just doing a GLM on the transformed covariates (if we ignore the penalization component). So are GAMs still linear in the same sense as GLMs? ($\eta = \sum \beta g(x)$?)*

### Qualitative Covariates / Factors

> This section is interesting because it shows how naturally linear models evolve into something like a mixed models through factor variables, and scratches the surface of how this all relates to ANOVA models.

*M&N pg. 52*  
Sets of observations are often indexed by one or more "classifying factors". Each factor is made up of a number of "indexes", or levels, which partition the data into disjoint groups or classes.  
For example, you may have data observing continuous variables wealth and age, and also have a education factor, which has levels "high school", "college", and "graduate".  
In a classic field experiment, a factor might be which crop variety is planted in a given plot.  

Factors may be of primary interest - e.g., if we have a "treatment" factor that partitions observations into treated and untreated groups, we are interested in the effect of treatment. Factors may also be secondary - included by necessity in the model, e.g., [blocking factors](https://en.wikipedia.org/wiki/Blocking_(statistics)).

Recall our linear predictor $\eta = \alpha + \beta x$.  
If $A$ is a factor with index $i$ (i.e., $i$ represents the levels of the factor), the model can be rewritten $\eta_i = \alpha_i + \beta x$
- Here the actual factor *becomes part of the intercept*, hence $\alpha_i$. This becomes obvious to me if we think of encoding the factor levels as dummies - say $A_i$ has 2 levels and we rewrite this as $\eta = \alpha + \gamma_1 A_1 + \gamma_2 A_2 + \beta x$. When an observation has $A_i = A_1$, the intercept becomes $\alpha + \gamma_1$. So we are shifting the intercept $\alpha$ based on the factor level.  

Note that this implies a separate intercept for each leavel of $A$, but a common slope $\beta$ assumed constant over the levels of the factor.

Often observations are "cross-classified" by many factors. If $A$, $B$, and $C$ are three such factors with indices $i$, $j$, and $k$, the simplest model ordinarily considered has form $\alpha_i + \beta_j + \gamma_k$.
- "This is the so-called main-effects model, which implies that if we arrange the data in a rectangular block and then look at cross-sections for each level of $A$, we'll find that it can be modelled by effects of $B$ and $C$ that are additive and equal in each cross-section... and similairly for the other factors."

- To achieve satisfactory fit, it may be necessary to include **interactions** between factors, e.g. $(\alpha \beta)_{ij}$ implies a separate effect for each combination of the indices $i$ and $j$ (a "two-factor interaction").

- If $i$ is the index for the levels of a factor $A$, we can represent the factor by a set of **dummy variates** so that the (intercept) term $\alpha_i$ can be written as $\alpha_1 u_1 + \dots + \alpha_k u_{k}$, where $k$ is the number of levels of the factor. (Terms "incidence vector" and "indicator vector" are also used.)  

### Mixed terms

*M&N pg. 55*  

Above we showed that our model $\eta_i = \alpha_i + \beta x$ has a separate intercept for each level of the factor $A$, but a constant slope $\beta$ over all levels.  
We may also consider a model where the slope is allowed to vary by level of the factor, requiring $\beta x$ to be replaced by $\beta_i x$.  
**Terms in the linear predictor that involve both continuous and qualitative variables are called mixed terms**.

M&N note that if one is to include a factor variable and not a mixed term, so assuming a constant slope, this assumption (of constant slope) should be tested - in simplest form, by comparing the fit of the model with constant slope with the fit when the slope is allowed to vary by level. (And therefore, stat software must provide a way to fit mixed terms).

Dummy variates for mixed terms take the same form as those for factors, but instead of 1s and 0s, the 1s are replaced with the value of the continuous variable.

## Model formulae

*M&N pg. 56*  

They describe the notation introduced by Wilkinson and Rogers (1973) for model formulae. It's interesting because it shows how the model formulae we use in R, etc are not just arbitrary syntax, but have a history and a logic behind them.

## Aliasing
We have a set of $n$-dimensional predictors $x_1, \dots, x_p$.  
These can be thought of as a set of $p$ vectors in $n$-dimensional Euclidean space. Then they define a subspace of *up to* $p$ dimensions - the maximum achieved only of the $x$'s are linearly independent. If $k$ linear relations exist among the $x$'s, the spanned subspace will have dimension $p - k$.



## Why we leave out a factor level when including factor covariates

*M&N pg. 64*  

When including factors in a regression, R defaults to dropping the first level, so that the coefficients (assuming Normal theory linear regression) on the levels of the factor are interpreted as measures of differences between each level's mean and the dropped level's mean. This is necessary to avoid "intrinsic aliasing" (aliasing due to intrinsic aspects of the model formula as opposed to idiosyncracies in the data):   
- Consider a model that contains an intercept and a factor variable (denoted as $1 + A$ in formula notation), where each level of the factor is represented as a dummy vector (so $y = \beta_0 + \beta_1 u_1 + \dots + \beta_k u_k$, for $k$ factor levels, where $u$ are dummy vectors where an element is 1 if the observation is in the factor level and 0 otherwise). An observation belongs to only one factor level, and so $\sum_{i=1}^k u_i = 1$ (i.e., if you sum across the rows of the dummy vectors, you have a vector of all 1s). This is a problem, because it means that the intercept (i.e., the vector of all 1s) is completely linearly dependent on the $u$'s, since it can be reproduced by simply adding them up (i.e., the intercept is intrinsically aliased).  
- To avoid this problem, we can put some "constraints" on the problem (constraints are not to be thought of as part of the model specification, but rather as a convenient way of resolving some ambiguity). Some (out of many) possible constraints given by M&N are:  
    - 1. Set the intercept $\beta_0$ to be zero so that each coefficient $\beta_1, \dots, \beta_k$ gives the of the corresponding group (in R, you could estimate a model without an intercept)
    - 2. Consider $\beta_1$ to be $0$, so that the intercept corresponds to the group mean of the first level and the other coefficients measure the differences between other group means and the first (this is the default in R)


## Why you should include the linear term if including polynomial terms, and the main effects if including interactions

*M&N pg. 69*

Even if covariates are not linearly related (i.e., they are not linearly dependent / "aliased") they may be functionally related. E.g., polynomial regression: $\beta_0 + \beta_1 x + \beta_2 x^2 + \dots$

Consider the formula $\beta_0 + \beta_1 x + \beta_2 x^2$.  
- If we instead fit $\beta_0 + \beta_2 x^2$, while excluding the linear term, this would imply that the maximum or minimum of the response occurs at $x = 0$, specifically (to see this, graph $y = cx^2$ and $y = bx + cx^2$ for different values of $b$ and $c$ - without the linear term, $y = cx^2$ has a min/max at 0, always)
- Since there is typically no reason to assume such a functional form, it doesn't make sense to exclude the linear term.  

Now consider the formula $\beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_3 x_1 x_2$, where $\beta_1 x_1$ and $\beta_2 x_2$ are considered main effects and $\beta_3 x_1 x_2$ is their interaction.  
- If we include the interaction term but not the main effects, this is equivalent to assuming that the point $(0,0)$ is a saddle-point of the response surface.
    - (again, simply graph this using a 3d graphing calculator. $z = xy$ is a surface with a saddle point at the origin. $z = x + xy$ moves the saddle point minus 1 in the y direction. $z = y + xy$ moves the saddle point minus 1 in the x direction. $z = x + y + xy$ has saddle point at $(-1, -1, -1)$. Adding coefficients other than 1 thus allows the saddle point to be moved anywhere, based on the data)
- Again, there is no reason to assume such a special property for the origin, so **the linear terms must be included with the cross-term**



## Selection of Covariates

*M&N pg. 89*  
- "Apart from the choice of link function and error distribution, the problem of modeling reduces to finding one or more appropriate parsimonious sets of covariates corresponding to a model matrix X..."
    - Parsimony implies that the model should be as simple as possible while still capturing the important features of the data. If a covariate has no detectable effect on the response, it should not be included
- At a minimum, the model should make sense physically - as state above, if interactions or polynomial terms are included, the main effects or lower-term polynomials should be included as well.
- Including a new covariate in a model is a balancing act. On one hand, including it may improve the fit - a reduction in the discrepancy between the data and fitted values. On the other hand, including it complicates the model (which may be undesirable for interpretability, *or in prediction sense may lead to overfitting... bias-variance tradeoff...*).



## On the asymptomtic distribution of the Deviance for the Binomial Distribution
*M&N pg. 118*

Recall that the residual deviance is defined to be twice the difference between the maximum achievable log-likelihood and the log-likelihood of the fitted model.  

The deviance function for the binomial distribution is:  
$D(Y; \hat{\pi}) = 2 [\sum_{i=1}^n y_i \log( \frac{y_i}{\hat{\mu}_i} ) + (m_i - y_i) \log ( \frac{m_i - y_i}{m_i - \hat{\mu}_i} )]$


It is often claimed that the random variable $D(Y; \hat{\pi})$ is asymptotically or approximately distributed as $\chi^2_{n-p}$, which is used to justify the use of $D$ as a goodness-of-fit statistic for testing the adequacy of the model.  
Proofs of the limiting $\chi^2_{n-p}$ distribution are based on the following specific assumptions:  
- The observations are distributed independently according to the binomial distribution - over-dispersion is not considered.  
- The approximation is based on a limiting operation in which the dimension of $Y$ is fixed at $n$ and $m_i \to \infty$ for each $i$ - and in fact, $m_i \pi_i (1 - \pi_i) \to \infty$.  
    - Under this limit, $D$ is approximately independent of the estimated parameters $\hat{\beta}$, and hence approximately independent of the fitted probabilities $\hat{\pi}$. Approximate independence is essential for $D$ to be considered as a g.o.f. statistic, but this property does not guarantee good power.  

If instead $n$ is large and $m_i \pi_i (1 - \pi_i)$ remains bounded, this theory breaks down:  
- The limiting $\chi^2_{n-p}$ approximation no longer holds. 
- $D$ is not independent of $\hat{\pi}$ even approximately, so a large value of $D$ could be obtained with high probability by judicious choice of $\beta$ and $\pi$.


The deviance function is most directly useful not as a direct measure of fit, but for comparing 2 nested models...
(discussed in notes 01).  
