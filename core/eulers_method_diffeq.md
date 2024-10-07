# Intro to Differential Equations: Solving separable equations, Euler’s
method, and direction fields
H Elliott

- [Euler’s method function](#eulers-method-function)
- [Example 1](#example-1)
- [Example 2](#example-2)
- [Some first-order ODEs and their direction fields from
  applications](#some-first-order-odes-and-their-direction-fields-from-applications)
  - [Newton’s law of cooling or
    heating](#newtons-law-of-cooling-or-heating)
  - [Infection in a population](#infection-in-a-population)
  - [Some random direction fields](#some-random-direction-fields)

Some basic examination of differential equations and their solutions.

Also a brief demonstration of Euler’s method for approximation of
solutions.

# Euler’s method function

``` r
#' Perform Euler's method to approximate the solution to an initial value problem
#' 
#' @param f The function defining the differential equation dy/dt = f(t, y), with
#'     signature `function(t, y)`
#' @param y0 The initial value of y at t0
#' @param t0 The initial value of t
#' @param t1 The final value of t
#' @param h The step size to take from t0 to t1
eulers <- function(f, y0, t0, t1, h) {
    t <- seq(t0, t1, by = h)
    y <- numeric(length(t))
    y[1] <- y0
    t[1] <- t0
    for (i in 2:length(t)) {
        y[i] <- y[i - 1] + h * f(t[i - 1], y[i - 1])
    }
    return(list(t = t, y = y, f = f, y0 = y0, t0 = t0, t1 = t1, h = h))
}
```

<details>
<summary>Code</summary>

``` r
# plotting functions
plot_eulers <- function(eulers_result, true_f = NULL, legend_pos = "topleft") {
    plot(eulers_result$t, eulers_result$y,
         type = "p", cex = 0.5, col = rgb(0, 0, 0, 0.5),
         main = "Euler's method results",
         sub = paste(deparse(eulers_result$f), collapse = ""),
         xlab = "t", ylab = "y")
    lines(eulers_result$t, eulers_result$y, col = rgb(0, 0, 0, 0.5))
    if (! is.null(true_f)) {
        tlin <- seq(res$t0, res$t1, by = 0.01)
        lines(tlin, true_f(tlin), col = "red")
    }
    if (legend_pos != "none") {
        legend(legend_pos, legend = c("Euler's method", "Exact solution"),
               col = c("black", "red"),
               lty = 1)
    }
}


#' plot direction field
#' f: differential equation function with signature `function(t, y)`
plot_directionfield <- function(f,
                                t_range = c(-5, 5),
                                y_range = c(-5, 5),
                                radius = 0.1, grid.by = 0.25, alpha = 1){
    # credit:
    # https://stackoverflow.com/questions/47984874/how-to-create-a-slope-field-in-r
    
    # initial plot - ensure large enough
    plot(t_range, y_range,
         main = "Direction field", ylab = "y", xlab = "t",
         pch = ".")
    # plot arrows
    tlin = seq(min(t_range), max(t_range), grid.by)
    ylin = seq(min(y_range), max(y_range), grid.by)
    for(x in tlin) {
        for(y in ylin) {
            slope = f(x, y)
            if(is.na(slope)) {
                col = rgb(0, 0, 0, alpha)
            } else if(slope > 0) {
                col = rgb(0, 0, 1, alpha)
            } else if (slope < 0) {
                col = rgb(1, 0, 0, alpha)
            } else if(slope == 0) {
                col = rgb(0, 1, 0, alpha)
            }
            arrows(radius * cos(atan(slope) + pi) + x,
                   radius * sin(atan(slope) + pi) + y,
                   radius * cos(atan(slope)) + x,
                   radius * sin(atan(slope)) + y, 
                   length = 0.2 * radius, col = col)
        }
    }
}
```

</details>

# Example 1

The differential equation  

![\frac{dy}{dt} = t (4 - y)](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%20t%20%284%20-%20y%29 "\frac{dy}{dt} = t (4 - y)")

has general solution

![y = 4 + C e^{-t^2/2}](https://latex.codecogs.com/svg.latex?y%20%3D%204%20%2B%20C%20e%5E%7B-t%5E2%2F2%7D "y = 4 + C e^{-t^2/2}")

where ![C](https://latex.codecogs.com/svg.latex?C "C") is a constant.

Derivation:

- ![\frac{dy}{4-y} = t dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7B4-y%7D%20%3D%20t%20dt "\frac{dy}{4-y} = t dt")  
- ![\int \frac{dy}{4-y} = \int t dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7Bdy%7D%7B4-y%7D%20%3D%20%5Cint%20t%20dt "\int \frac{dy}{4-y} = \int t dt")  
- ![- \ln\|4 - y\| = \frac{t^2}{2} + C_0](https://latex.codecogs.com/svg.latex?-%20%5Cln%7C4%20-%20y%7C%20%3D%20%5Cfrac%7Bt%5E2%7D%7B2%7D%20%2B%20C_0 "- \ln|4 - y| = \frac{t^2}{2} + C_0")  
- ![\ln\|4 - y\| = -\frac{t^2}{2} + C_1](https://latex.codecogs.com/svg.latex?%5Cln%7C4%20-%20y%7C%20%3D%20-%5Cfrac%7Bt%5E2%7D%7B2%7D%20%2B%20C_1 "\ln|4 - y| = -\frac{t^2}{2} + C_1")  
- ![\|4 - y\| = e^{-t^2/2} e^{C_1}](https://latex.codecogs.com/svg.latex?%7C4%20-%20y%7C%20%3D%20e%5E%7B-t%5E2%2F2%7D%20e%5E%7BC_1%7D "|4 - y| = e^{-t^2/2} e^{C_1}")  
- ![4 - y = \pm e^{-t^2/2} e^{C_1}](https://latex.codecogs.com/svg.latex?4%20-%20y%20%3D%20%5Cpm%20e%5E%7B-t%5E2%2F2%7D%20e%5E%7BC_1%7D "4 - y = \pm e^{-t^2/2} e^{C_1}")  
- ![y = 4 \pm e^{-t^2/2} e^{C_1}](https://latex.codecogs.com/svg.latex?y%20%3D%204%20%5Cpm%20e%5E%7B-t%5E2%2F2%7D%20e%5E%7BC_1%7D "y = 4 \pm e^{-t^2/2} e^{C_1}")  
- ![y = 4 + C e^{-t^2/2}](https://latex.codecogs.com/svg.latex?y%20%3D%204%20%2B%20C%20e%5E%7B-t%5E2%2F2%7D "y = 4 + C e^{-t^2/2}")
  where
  ![C = \pm e^{C_1}](https://latex.codecogs.com/svg.latex?C%20%3D%20%5Cpm%20e%5E%7BC_1%7D "C = \pm e^{C_1}")
  is a constant.

One solution to the differential equation is
![y = 4](https://latex.codecogs.com/svg.latex?y%20%3D%204 "y = 4") so
that
![\frac{dy}{dt} = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%200 "\frac{dy}{dt} = 0")
for all ![t](https://latex.codecogs.com/svg.latex?t "t"). This is
satisfied if
![C = 0](https://latex.codecogs.com/svg.latex?C%20%3D%200 "C = 0").
Other particular solutions depend on the initial condition.

If
![y(0) = 0](https://latex.codecogs.com/svg.latex?y%280%29%20%3D%200 "y(0) = 0"),
then
![C = -4](https://latex.codecogs.com/svg.latex?C%20%3D%20-4 "C = -4")
and the solution is

![y = 4 - 4 e^{-t^2/2}](https://latex.codecogs.com/svg.latex?y%20%3D%204%20-%204%20e%5E%7B-t%5E2%2F2%7D "y = 4 - 4 e^{-t^2/2}")

``` r
f <- \(t, y) t * (4 - y)
true_sol <- \(t) 4 - 4 * exp( (-t^2) / 2) # if y(t = 0) = 0

par(mfrow = c(2, 2))
for (stepsize in c(0.5, 0.1, 0.05, 0.01)) {
    res <- eulers(f, y0 = 0, t0 = 0, t1 = 4, h = stepsize)
    plot_eulers(res, true_sol, legend_pos = "none")
    abline(h = 4, col = "blue", lty = 2)
    mtext(side=3, line=0.5, at=-0.01, adj=0, cex=0.7,
          paste("step size = ", res$h))
}
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-3-1.png)

Note that the red lines show the exact solution and the black
points/lines show the approximated solution via Euler’s method.

A different initial condition:

If
![y(0) = 8](https://latex.codecogs.com/svg.latex?y%280%29%20%3D%208 "y(0) = 8"),
then ![C = 4](https://latex.codecogs.com/svg.latex?C%20%3D%204 "C = 4")
and the solution is

![y = 4 + 4 e^{-t^2/2}](https://latex.codecogs.com/svg.latex?y%20%3D%204%20%2B%204%20e%5E%7B-t%5E2%2F2%7D "y = 4 + 4 e^{-t^2/2}")

``` r
f <- \(t, y) t * (4 - y)
true_sol <- \(t) 4 + 4 * exp( (-t^2) / 2) # if y < 4 and y(t = 0) = 0

par(mfrow = c(2, 2))
for (stepsize in c(0.5, 0.1, 0.05, 0.01)) {
    res <- eulers(f, y0 = 8, t0 = 0, t1 = 4, h = stepsize)
    plot_eulers(res, true_sol, legend_pos = "none")
    abline(h = 4, col = "blue", lty = 2)
    mtext(side=3, line=0.5, at=-0.01, adj=0, cex=0.7,
          paste("step size = ", res$h))
}
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-4-1.png)

Notice that in both cases, the solutions seems to converge to
![y = 4](https://latex.codecogs.com/svg.latex?y%20%3D%204 "y = 4") (the
blue dashed line) as ![t](https://latex.codecogs.com/svg.latex?t "t")
increases. This is because
![y = 4](https://latex.codecogs.com/svg.latex?y%20%3D%204 "y = 4") is a
stable solution to the differential equation.  
This can be seen by plotting the direction field of the differential
equation.  
The black lines plot the particular solutions examined above. Each
solution is tangent to the direction field at every point.

``` r
plot_directionfield(f,
                    t_range = c(-5, 5), y_range = c(0, 8),
                    alpha = 0.5)
mtext(side=3, line=0.5, at=-3, adj=0, cex=0.7,
      "dy/dt = t(4 - y), and solutions for y(0) = 0 and y(0) = 8")

tlin <- seq(-5, 5, 0.01)
lines(tlin, 4 - 4 * exp( (-tlin^2) / 2) , col = "black", lty = 1)
lines(tlin, 4 + 4 * exp( (-tlin^2) / 2) , col = "black", lty = 1)
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-5-1.png)

# Example 2

Here is another simple example.  
Differential equation

![\frac{dy}{dt} = y](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%20y "\frac{dy}{dt} = y")

has general solution

![y = C e^t](https://latex.codecogs.com/svg.latex?y%20%3D%20C%20e%5Et "y = C e^t")

where ![C](https://latex.codecogs.com/svg.latex?C "C") is a constant:

- ![\frac{dy}{y} = dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7By%7D%20%3D%20dt "\frac{dy}{y} = dt")  
- ![\int \frac{dy}{y} = \int dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7Bdy%7D%7By%7D%20%3D%20%5Cint%20dt "\int \frac{dy}{y} = \int dt")  
- ![\ln\|y\| = t + C_0](https://latex.codecogs.com/svg.latex?%5Cln%7Cy%7C%20%3D%20t%20%2B%20C_0 "\ln|y| = t + C_0")  
- ![\|y\| = e^{t + C_0}](https://latex.codecogs.com/svg.latex?%7Cy%7C%20%3D%20e%5E%7Bt%20%2B%20C_0%7D "|y| = e^{t + C_0}")
- ![y = \pm e^{t}e^{C_0}](https://latex.codecogs.com/svg.latex?y%20%3D%20%5Cpm%20e%5E%7Bt%7De%5E%7BC_0%7D "y = \pm e^{t}e^{C_0}")  
- ![y = C e^t](https://latex.codecogs.com/svg.latex?y%20%3D%20C%20e%5Et "y = C e^t")
  where
  ![C = \pm e^{C_0}](https://latex.codecogs.com/svg.latex?C%20%3D%20%5Cpm%20e%5E%7BC_0%7D "C = \pm e^{C_0}")
  is a constant.

Notice that if If
![y = 0](https://latex.codecogs.com/svg.latex?y%20%3D%200 "y = 0") then
![\frac{dy}{dt} = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%200 "\frac{dy}{dt} = 0")
for all ![t](https://latex.codecogs.com/svg.latex?t "t"), and this is
achieved when
![C = 0](https://latex.codecogs.com/svg.latex?C%20%3D%200 "C = 0").
Otherwise the particular solution depends on the initial condition.

If
![y(t = 0) = 2](https://latex.codecogs.com/svg.latex?y%28t%20%3D%200%29%20%3D%202 "y(t = 0) = 2"),
![C = 2](https://latex.codecogs.com/svg.latex?C%20%3D%202 "C = 2") and
the particular solution is
![y = 2 e^t](https://latex.codecogs.com/svg.latex?y%20%3D%202%20e%5Et "y = 2 e^t").  
Generally, for this simple equation if we start at
![t = 0](https://latex.codecogs.com/svg.latex?t%20%3D%200 "t = 0"),
whatever the value of
![y(0)](https://latex.codecogs.com/svg.latex?y%280%29 "y(0)") is will be
the value of ![C](https://latex.codecogs.com/svg.latex?C "C"). Here are
some plots for different initial conditions.

``` r
f <- \(t, y) y
true_sol <- \(t, y0) y0 * exp(t)

par(mfrow = c(2, 2))
tlin <- seq(0, 4, by = 0.01)
for (y0 in c(-2, -1, 1, 2)) {
    plot(tlin, true_sol(tlin, y0), type = "l", col = "red",
         main = paste("Exact solution when y(0) = ", y0),
         sub = paste("y =", y0, "* exp(t)"),
         xlab = "t", ylab = "y")
    abline(h = 0, col = "blue", lty = 2)
}
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-6-1.png)

In this case, the solution
![y = 0](https://latex.codecogs.com/svg.latex?y%20%3D%200 "y = 0") is
unstable because solutions that start near
![y = 0](https://latex.codecogs.com/svg.latex?y%20%3D%200 "y = 0") will
tend move away from it as
![t](https://latex.codecogs.com/svg.latex?t "t") increases.

``` r
plot_directionfield(f,
                    t_range = c(-5, 5), y_range = c(-5, 5),
                    alpha = 0.5)
mtext(side=3, line=0.5, at=-3, adj=0, cex=0.7,
      "dy/dt = y, and solutions for y(0) = 1 and y(0) = -1")

tlin <- seq(-5, 5, 0.01)
lines(tlin, 1 * exp(tlin) , col = "black", lty = 1)
lines(tlin, -1 * exp(tlin) , col = "black", lty = 1)
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-7-1.png)

# Some first-order ODEs and their direction fields from applications

## Newton’s law of cooling or heating

The change in the temperature of an object, which we will denote as
![y](https://latex.codecogs.com/svg.latex?y "y"), over time,
![t](https://latex.codecogs.com/svg.latex?t "t") converges to the
temperature of its surroundings, or “room temperature”,
![R](https://latex.codecogs.com/svg.latex?R "R").  
The parameter ![k](https://latex.codecogs.com/svg.latex?k "k") is the
“cooling rate”, which depends on the object.

The differential equation is:

![\frac{dy}{dt} = -k (y - R)](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%20-k%20%28y%20-%20R%29 "\frac{dy}{dt} = -k (y - R)")

Solving:

- ![\frac{-1}{k} \frac{dy}{R - y} = dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7B-1%7D%7Bk%7D%20%5Cfrac%7Bdy%7D%7BR%20-%20y%7D%20%3D%20dt "\frac{-1}{k} \frac{dy}{R - y} = dt")  
- ![\int \frac{-1}{k} \frac{dy}{R - y} = \int dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7B-1%7D%7Bk%7D%20%5Cfrac%7Bdy%7D%7BR%20-%20y%7D%20%3D%20%5Cint%20dt "\int \frac{-1}{k} \frac{dy}{R - y} = \int dt")  
- ![\frac{-1}{k} \ln\|R - y\| = t + C_0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B-1%7D%7Bk%7D%20%5Cln%7CR%20-%20y%7C%20%3D%20t%20%2B%20C_0 "\frac{-1}{k} \ln|R - y| = t + C_0")  
- ![\ln\|R - y\| = -k t + C_1](https://latex.codecogs.com/svg.latex?%5Cln%7CR%20-%20y%7C%20%3D%20-k%20t%20%2B%20C_1 "\ln|R - y| = -k t + C_1")  
- ![\|R - y\| = e^{-kt + C_1}](https://latex.codecogs.com/svg.latex?%7CR%20-%20y%7C%20%3D%20e%5E%7B-kt%20%2B%20C_1%7D "|R - y| = e^{-kt + C_1}")  
- ![R - y = \pm e^{C_1} e^{-kt}](https://latex.codecogs.com/svg.latex?R%20-%20y%20%3D%20%5Cpm%20e%5E%7BC_1%7D%20e%5E%7B-kt%7D "R - y = \pm e^{C_1} e^{-kt}")  
- ![y = R \pm e^{C_1} e^{-kt}](https://latex.codecogs.com/svg.latex?y%20%3D%20R%20%5Cpm%20e%5E%7BC_1%7D%20e%5E%7B-kt%7D "y = R \pm e^{C_1} e^{-kt}")  
- ![y = R + A e^{-kt}](https://latex.codecogs.com/svg.latex?y%20%3D%20R%20%2B%20A%20e%5E%7B-kt%7D "y = R + A e^{-kt}")
  where
  ![A = \pm e^{C_1}](https://latex.codecogs.com/svg.latex?A%20%3D%20%5Cpm%20e%5E%7BC_1%7D "A = \pm e^{C_1}")
  is a constant.

Generally, if the initial condition is that
![y(t = 0) = y_0](https://latex.codecogs.com/svg.latex?y%28t%20%3D%200%29%20%3D%20y_0 "y(t = 0) = y_0"),
then
![y_0 = R + Ae^0](https://latex.codecogs.com/svg.latex?y_0%20%3D%20R%20%2B%20Ae%5E0 "y_0 = R + Ae^0")
and so
![A = y_0 - R](https://latex.codecogs.com/svg.latex?A%20%3D%20y_0%20-%20R "A = y_0 - R").

So if the initial temperature of the object,
![y_0](https://latex.codecogs.com/svg.latex?y_0 "y_0"), is greater than
the room temperature, the object will cool down. Otherwise, it will heat
up. And when
![y = R](https://latex.codecogs.com/svg.latex?y%20%3D%20R "y = R"), the
object has reached the room temperature, and from the differential
equation we can see that
![\frac{dy}{dt} = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%200 "\frac{dy}{dt} = 0"),
so the object will stay at room temp.

``` r
# Newton's law of cooling or heating
# y = temperature of the object
# R = room temperature/temperature of the surroundings
# k = the cooling rate (object-dependent), k > 0
make_diffeq <- \(R, k) return( \(t, y) -k * (y - R) )

# plot direction field
plot_directionfield(f = make_diffeq(R = 5, k = 0.5),
                    t_range = c(0, 10), y_range = c(4, 8))


sol <- \(t, A, R, k) R + A * exp(-k * t)
# y(0) = 10, R = 5 # (object is cooling down to R)
tlin <- seq(0, 10, 0.01)
lines(tlin, sol(tlin, A = 5, R = 5, k = 0.5), col = "black", lty = 1)
# y(0) = 0, R = 5  # (object is warming up to R)
tlin <- seq(0, 10, 0.01)
lines(tlin, sol(tlin, A = -5, R = 5, k = 0.5), col = "black", lty = 1)
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-8-1.png)

## Infection in a population

The change in the size of the infected population,
![y](https://latex.codecogs.com/svg.latex?y "y"), over time,
![t](https://latex.codecogs.com/svg.latex?t "t"), with total population,
![P](https://latex.codecogs.com/svg.latex?P "P"), is given by:

![\frac{dy}{dt} = k y (P - y)](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%20k%20y%20%28P%20-%20y%29 "\frac{dy}{dt} = k y (P - y)")

This can be though of as being proportional to the number of infected
people (![y](https://latex.codecogs.com/svg.latex?y "y")) times the
number of uninfected people in the population
(![P - y](https://latex.codecogs.com/svg.latex?P%20-%20y "P - y")),
where ![k](https://latex.codecogs.com/svg.latex?k "k") is a constant
that acts like an “infection rate” for the particular disease.

This is a separable equation so it can be solved - you need to use
partial fractions to integrate (steps shown below).  
The solution:

![y = \frac{P C e^{Pkt}}{1 + C e^{Pkt}}](https://latex.codecogs.com/svg.latex?y%20%3D%20%5Cfrac%7BP%20C%20e%5E%7BPkt%7D%7D%7B1%20%2B%20C%20e%5E%7BPkt%7D%7D "y = \frac{P C e^{Pkt}}{1 + C e^{Pkt}}")

The particular solutions are logistic curves, as you can tell by their
S-shape.

![C](https://latex.codecogs.com/svg.latex?C "C") is a constant that will
depend on the initial condition.

If
![y(t = 0) = 1](https://latex.codecogs.com/svg.latex?y%28t%20%3D%200%29%20%3D%201 "y(t = 0) = 1")
(one infected person at time 0), then
![1 = \frac{P C}{1 + C}](https://latex.codecogs.com/svg.latex?1%20%3D%20%5Cfrac%7BP%20C%7D%7B1%20%2B%20C%7D "1 = \frac{P C}{1 + C}")
so
![C = \frac{1}{P - 1}](https://latex.codecogs.com/svg.latex?C%20%3D%20%5Cfrac%7B1%7D%7BP%20-%201%7D "C = \frac{1}{P - 1}").

``` r
make_diffeq <- \(P, k) return( \(t, y) k * y * (P - y) )

# solution when y(0) = 1:
sol <- \(t, P, k, C) {
    e_pkt <- exp(P * k * t)
    (P * C * e_pkt) / (1 + C * e_pkt)
} 

# plot direction field
plot_directionfield(f = make_diffeq(P = 10, k = 0.1),
                    t_range = c(0, 10), y_range = c(0, 12))
mtext(side=3, line=0.5, at=0, adj=0, cex=0.7,
      paste("Initial population =", 10))

# plot solution for a few different infection rates
tlin <- seq(0, 10, 0.01)
lines(tlin, sol(tlin, P = 10, k = 0.05, C = 1/9), col = "black", lty = 1)
lines(tlin, sol(tlin, P = 10, k = 0.1, C = 1/9), col = "black", lty = 1)
lines(tlin, sol(tlin, P = 10, k = 0.3, C = 1/9), col = "black", lty = 1)
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-9-1.png)

Solutions for three different infection rates
(![k](https://latex.codecogs.com/svg.latex?k "k")) are shown - the
larger the ![k](https://latex.codecogs.com/svg.latex?k "k"), the faster
the infected population ![y](https://latex.codecogs.com/svg.latex?y "y")
grows to the total population
![P](https://latex.codecogs.com/svg.latex?P "P").

Solving the differential equation:

- ![\frac{dy}{dt} = k y (P - y)](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7Bdt%7D%20%3D%20k%20y%20%28P%20-%20y%29 "\frac{dy}{dt} = k y (P - y)")  
- ![\frac{dy}{y (P - y)} = k dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdy%7D%7By%20%28P%20-%20y%29%7D%20%3D%20k%20dt "\frac{dy}{y (P - y)} = k dt")  
- ![\int \frac{dy}{y (P - y)} = \int k dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7Bdy%7D%7By%20%28P%20-%20y%29%7D%20%3D%20%5Cint%20k%20dt "\int \frac{dy}{y (P - y)} = \int k dt")
  - ![\frac{1}{y(P - y)} = \frac{A(P -y) + By}{y(P-y)} = \frac{(B-A)y + AP}{y(P-y)}](https://latex.codecogs.com/svg.latex?%5Cfrac%7B1%7D%7By%28P%20-%20y%29%7D%20%3D%20%5Cfrac%7BA%28P%20-y%29%20%2B%20By%7D%7By%28P-y%29%7D%20%3D%20%5Cfrac%7B%28B-A%29y%20%2B%20AP%7D%7By%28P-y%29%7D "\frac{1}{y(P - y)} = \frac{A(P -y) + By}{y(P-y)} = \frac{(B-A)y + AP}{y(P-y)}")  
  - so
    ![B - A = 0 \Rightarrow A = B](https://latex.codecogs.com/svg.latex?B%20-%20A%20%3D%200%20%5CRightarrow%20A%20%3D%20B "B - A = 0 \Rightarrow A = B")  
  - and
    ![AP = 1 \Rightarrow A = B = \frac{1}{P}](https://latex.codecogs.com/svg.latex?AP%20%3D%201%20%5CRightarrow%20A%20%3D%20B%20%3D%20%5Cfrac%7B1%7D%7BP%7D "AP = 1 \Rightarrow A = B = \frac{1}{P}")  
  - so
    ![\int \frac{dy}{y (P-y)} = 1/P \int dy/y + 1/P \int dy/(P-y) dy](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7Bdy%7D%7By%20%28P-y%29%7D%20%3D%201%2FP%20%5Cint%20dy%2Fy%20%2B%201%2FP%20%5Cint%20dy%2F%28P-y%29%20dy "\int \frac{dy}{y (P-y)} = 1/P \int dy/y + 1/P \int dy/(P-y) dy")  
- ![\frac{1}{P} \ln\|y/(P-y)\| = kt + C_0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B1%7D%7BP%7D%20%5Cln%7Cy%2F%28P-y%29%7C%20%3D%20kt%20%2B%20C_0 "\frac{1}{P} \ln|y/(P-y)| = kt + C_0")  
- ![\ln\|y/(P-y)\| = Pkt + C_0P](https://latex.codecogs.com/svg.latex?%5Cln%7Cy%2F%28P-y%29%7C%20%3D%20Pkt%20%2B%20C_0P "\ln|y/(P-y)| = Pkt + C_0P")  
- ![y/(P-y) = \pm e^{Pkt + C_0P} = y / (P - y) = Ce^{Pkt}](https://latex.codecogs.com/svg.latex?y%2F%28P-y%29%20%3D%20%5Cpm%20e%5E%7BPkt%20%2B%20C_0P%7D%20%3D%20y%20%2F%20%28P%20-%20y%29%20%3D%20Ce%5E%7BPkt%7D "y/(P-y) = \pm e^{Pkt + C_0P} = y / (P - y) = Ce^{Pkt}"),
  where
  ![C = \pm e^{C_0P}](https://latex.codecogs.com/svg.latex?C%20%3D%20%5Cpm%20e%5E%7BC_0P%7D "C = \pm e^{C_0P}")  
- ![y = (P-y)Ce^{Pkt} = PC{e^{Pkt}} - yCe^{Pkt}](https://latex.codecogs.com/svg.latex?y%20%3D%20%28P-y%29Ce%5E%7BPkt%7D%20%3D%20PC%7Be%5E%7BPkt%7D%7D%20-%20yCe%5E%7BPkt%7D "y = (P-y)Ce^{Pkt} = PC{e^{Pkt}} - yCe^{Pkt}")  
- ![y + yCe^{Pkt} = PC{e^{Pkt}}](https://latex.codecogs.com/svg.latex?y%20%2B%20yCe%5E%7BPkt%7D%20%3D%20PC%7Be%5E%7BPkt%7D%7D "y + yCe^{Pkt} = PC{e^{Pkt}}")  
- ![y(1 + Ce^{Pkt}) = PC{e^{Pkt}}](https://latex.codecogs.com/svg.latex?y%281%20%2B%20Ce%5E%7BPkt%7D%29%20%3D%20PC%7Be%5E%7BPkt%7D%7D "y(1 + Ce^{Pkt}) = PC{e^{Pkt}}")  
- ![y = \frac{P C e^{Pkt}}{1 + C e^{Pkt}}](https://latex.codecogs.com/svg.latex?y%20%3D%20%5Cfrac%7BP%20C%20e%5E%7BPkt%7D%7D%7B1%20%2B%20C%20e%5E%7BPkt%7D%7D "y = \frac{P C e^{Pkt}}{1 + C e^{Pkt}}")

## Some random direction fields

``` r
par(mfrow = c(1, 2))
plot_directionfield(f = \(y, t) t + y)
plot_directionfield(f = \(t, y) sin(t * y))
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
par(mfrow = c(1, 2))
plot_directionfield(f = \(t, y) sin(t) * y)
plot_directionfield(f = \(t, y) sin(t) * y + cos(t))
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
par(mfrow = c(1, 2))
plot_directionfield(f = \(y, t) sin(t) + cos(y))
plot_directionfield(f = \(y, t) sin(t) + cos(y) + y)
```

![](eulers_method_diffeq_files/figure-commonmark/unnamed-chunk-12-1.png)
