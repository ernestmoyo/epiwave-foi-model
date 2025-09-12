library(greta)
library(DiagrammeR)
# data
x <- as_data(iris$Petal.Length)
y <- as_data(iris$Sepal.Length)

# variables

int <- normal(0,1)
coef <- normal(0,3)
sd <- student(3, 0, 1, truncation = c(0, Inf))

# operations
mean <- int + coef*x
#likelihood
distribution(y) <- normal(mean, sd)
#defining the model
m <- model(int, coef, sd)
#plotting
plot(m)
