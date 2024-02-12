
# Generalized Linear Models, continued

Notes based on [the book by Nelder and McCullagh, 1989](https://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf).


Some random-ish bits from the book that were interesting but not so directly related to deciphering GLMs.  


## Functions of Continuous Covariates

*M&N pg. 51*  
- In the context of GLMs, linearity refers to linearity of $\eta$ in terms of its parameters.  
- A continuous covariate $x$ in a model may be replaced by an arbitrary function $g(x)$ without destroying the linearity of the model. For example, we may add $x^2$ and $x^3$
    - *When reading about GAMs: after the transformation of a covariate into basis functions, we are essentially just doing a GLM on the transformed covariates. So are GAMs still linear in the same sense as GLMs? ($\eta = \sum \beta g(x)$?)*


## Why we leave out a factor level when including factor covariates

*M&N pg. 64*  
There is an interesting section that discusses model formulae and the different types of possible covariate combinations, with a note on aliasing - i.e., linear dependence or perfect multicollinearity.  

When including factors in a regression, R defaults to dropping the first level, so that the coefficients (assuming Normal theory linear regression) on the levels of the factor are interpreted as measures of differences between each level's mean and the dropped level's mean. This is necessary to avoid "intrinsic aliasing" (aliasing due to intrinsic aspects of the model formula as opposed to idiosyncracies in the data):   
- Consider a model that contains an intercept and a factor variable (denoted as $1 + A$ in their notation), where each level of the factor is represented as a dummy vector (so $y = \beta_0 + \beta_1 u_1 + \dots + \beta_k u_k$, for $k$ factor levels, where $u$ are dummy vectors where an element is 1 if the observation is in the factor level and 0 otherwise). An observation belongs to only one factor level, and so $\sum_{i=1}^k u_i = 1$ (i.e., if you sum across the rows of the dummy vectors, you have a vector of all 1s). This is a problem, because it means that the intercept (i.e., the vector of all 1s) is completely linearly dependent on the $u$'s, since it can be reproduced by simply adding them up (i.e., the intercept is intrinsically aliased).  
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