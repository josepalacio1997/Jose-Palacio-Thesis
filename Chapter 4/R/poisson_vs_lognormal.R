par(mfrow = c(2,2))

lambdas <- c(10, 20, 30, 50)

for(lambda in lambdas){
  
  mu_logn <- log(lambda) - 0.5 / lambda
  sigma_logn <- 1 / sqrt(lambda)
  
  x <- 0:(lambda + 5 * sqrt(lambda))
  
  # Densidad LogNormal
  plot(
    x,
    dlnorm(x, meanlog = mu_logn, sdlog = sigma_logn),
    type = "l",
    lwd = 2,
    col = "red",
    main = bquote(lambda == .(lambda)),
    xlab = "x",
    ylab = "Density / Probability"
  )
  
  # PMF Poisson
  points(
    x,
    dpois(x, lambda),
    pch = 16,
    col = "black"
  )
  
  legend(
    "topright",
    legend = c("LogNormal Approx", "Poisson"),
    col = c("red", "black"),
    lty = c(1, NA),
    pch = c(NA, 16),
    bty = "n"
  )
}