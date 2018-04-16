
Fm <- function(x,a,b,c,d,e,f,g) {
    x = log(x * 29)  ## transform to log years
    return (exp((a*x*x + b*x + c)*(d*x + e - cos( (x - f) / g ) )))
}

## model for ceu population size
## x = time in generations
ceu <- function(x) {
    x = pmin(4100000/29, x)
    a = 0.02869
    b = -0.5555
    c = 3.016
    d = -3.651
    e = 63.3
    f = 10.8
    g = 0.51
    return (Fm(x,a,b,c,d,e,f,g))
}

yri <- function(x) {
    return (ceu(pmax(5555,x)))
}

migr_ceu <- function(x) {
        ## version of 10/4/2018 (model 2): push back migration into ceu a little.  Looks good now.
    if (x < 3448) return (0)         # 100 kya
    if (x > 6900) return (0.00025)   # 200 kya
    return (0.00025 * (1.0 - (log(x)-log(6900))/(log(3448)-log(6900))))
}

## model for migration from ceu into yri, forward in time
migr_yri <- function(x, variant) {
    start = 1379
    mid = 2069
    end = 3103   ## 40k, 60k, 90k; "model 3"
    strength = 0.001
    if (variant == -1) {
        start = 1025
        mid = 1700
        end = 2900  ## # 30k, 50k, 85k; "model 1", 9/4/2018
    }
    if (variant == 1) {
        start = 1724
        mid = 2586
        end = 3448 ## 50k, 75k, 100k; "model 2", 10/4/2018, moved peak to around 75kya
    }
    if (x < start) return (0)
    if (x > end) return (0)
    if (x < mid) {
        return (strength * (1.0 - (log(x)-log(mid)) / (log(start)-log(mid))))
    }
    return (strength * (1.0 - (log(x)-log(mid)) / (log(end)-log(mid))))
}




    
