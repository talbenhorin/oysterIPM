library(tidyverse)
library(ggplot2)
out %>%
  ggplot(aes(x=time, y=I))+
  geom_line() +
  labs(title = "Fraction of population infected over time",
       subtitle = "50 year simulation - steady state is not reached",
       x = "Time",
       y = "Fraction Infected")

out%>%
  ggplot(aes(x = S, y=I))+
  geom_point(size = 0.5)+
  labs(x = "Percent Susceptible",
       y = "Percent Infected",
       title = "Joint Time Series")