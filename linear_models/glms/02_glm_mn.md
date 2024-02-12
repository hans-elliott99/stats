
# Generalized Linear Models, continued

Notes based on [the book by Nelder and McCullagh, 1989](https://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf).


# Normal/Gaussian GLMs with constant variance

When our outcome is continuous and its possible values span the whole real line, we can use the Normal distribution to model $Y|X$.  
In particular, if we assume a constant variance across observations, we can build a model like:  
- Random Component: $Y |X \sim N_n(\mu, \sigma^2 I_n)$
- Systematic Component: $\eta = X \beta$
- Identity Link: $\mu = \eta$  

Thus, $E[Y|X] = X \beta$.  
Thus we can make the same assumptions we make under classical ("Normal theory") linear regression and it can be shown that the classical linear regression is this case of the Normal GLM with constant variance.  
- Recall, in linear regression, minimizing least squares can be shown to be identical to likelihood maximization 

Arguably, the most important assumption here is that of constant variance across observations, since "the theory of least squares can be developed using only first- and second-moment assumptions in addition to independence, without requiring the additional assumption of Normality"
- This should be checked, for example graphically by examining the residuals

## When to use it

Obviously we can/should use this case of the GLM if we have an outcome variable that matches the assumptions of a Normal r.v. - e.g., if the outcome is continuous and can take any value on the real line.  

M&N expand this a bit though, saying:
- It can also be used as an approximation for discrete measurements
- It can be used to model continuous data that is strictly postiive (e.g., height or weight) as long as the data are sufficiently far from zero. For example if data have a mean of 100 and s.d. of 10, the part of the Normal distribution covering the negative half of the real line is is negligible for most practical purposes
    - If the data are or could be close to zero, then this would be a bad idea, and the gamma distribution might be used instead
- *I thought that this was interesting since GLMs are often mentioned as a required alternative to classic linear regression specifically to model data that isn't continuous or spanning all of R. It gives some more strength to the approach, in econometrics for example, where linear regression is used almost exclusively for causal inference-type modeling since it is relatively interpretable. If data is sufficiently far from 0, it very well may be worth modeling it as Normal, as opposed to something else. All models are wrong after all.*     

