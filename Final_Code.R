
library(dplyr)
library(ggplot2)
library(forecast)
library(lubridate)

# Reading the data
data1 <- read.csv("air-quality-india.csv")

# Making hourly observations into weekly average
data2 <- data1 %>%
  mutate(week = week(Timestamp)) %>%
  group_by(Year,week) %>%
  summarise_at(vars("PM2.5"), mean)

# Selecting data that starts from 1st week of 2018
data2<-data2[-c(1:9),]

#Building weekly time-series
weekly <- ts(data2$PM2.5, start = c(2018,1), end = c(2022,23), freq = 52)

# Observing decomposed model
plot(weekly)
plot(decompose(weekly))

    # It seems like there is mostly seasonality, very little trend 
    # The seasonality seems additive rather than multiplicative 
    # meaning it doesn't become larger of smaller with time.
    # Therefore, visually we should think about a seasonal naive model, 
    # or a linear model with seasonality, and maybe add some additive trend 
    # Nevertheless we will try multiple models to make sure

# Dividing Data set into train and validation
nValid <- 52          # There are 52 weeks in a year
nTrain <- length(weekly)-nValid          

train.ts <- window(weekly, start = c(2018,1), end = c(2018,nTrain))
valid.ts <- window(weekly, start = c(2018,nTrain+1), end = c(2018, nTrain+nValid))

# Observing validation and training data set plots (not necessary)
plot(train.ts)
plot(valid.ts)

#--------------------------------------------------------------------------------

# Building simple moving average model (not necessary since too basic model for this data)
moving_average <- ma(train.ts, order = 4)
ma_forecast <- forecast(moving_average, h = nValid)
plot(ma_forecast)
lines(valid.ts, col = "grey20", lty = 3) 

# Building simple simple mean model (not necessary since large seasonality present)
average_value <- meanf(train.ts, h = nValid)     
autoplot(average_value, facets = FALSE) + 
  xlab("Week") + ylab("Demand") + ggtitle("Time Series for the weekly demand") 
lines(valid.ts, col = "grey20", lty = 3) 

# Building Seasonal Naive models
seasonnal_naive <- snaive(train.ts, h = nValid)  
seas_naive_forecast <- forecast(seasonnal_naive, h = nValid)
accuracy(seas_naive_forecast, valid.ts)

# Multiple plots of Seasonal Naive since it seems like a pretty good model
# Plot of only the forecasted part (with CI) with the actual data in the background 
plot(seas_naive_forecast, xlim = c(2021.5,2022.5), ylab = "PM 2.5", xlab = "Years")  
lines(valid.ts, col = "grey20", lty = 3) 

# Plot all the forecasted part (with CI) with the actual data in the background 
plot(seas_naive_forecast, ylab = "PM 2.5", xlab = "Years")  
lines(valid.ts, col = "grey20", lty = 3) 
lines(seas_naive_forecast$fitted, col = "blue", lwd = 2)

# Plot the residuals of the forecast
plot(seas_naive_forecast$residuals, ylab = "Forecast Errors", xlab = "Years")

# Building linear model with season --> the one we believe it will be 
train.lm1 <- tslm(train.ts ~ season)
train.lm1.pred <- forecast(train.lm1, h = nValid)
accuracy(train.lm1.pred, valid.ts)
plot(train.lm1.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

# Building linear model with trend
train.lm2 <- tslm(train.ts ~ trend)
train.lm2.pred <- forecast(train.lm2, h = nValid)
accuracy(train.lm2.pred, valid.ts)
plot(train.lm2.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

# Building linear model with additive seasonality and trend
train.lm3 <- tslm(train.ts ~ season + trend)
train.lm3.pred <- forecast(train.lm3, h = nValid)
accuracy(train.lm3.pred, valid.ts)
plot(train.lm3.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

# Building linear model with additive seasonality and exponential trend
train.lm6 <- tslm(train.ts ~ season + (trend + I(trend^2)))
train.lm6.pred <- forecast(train.lm6, h = nValid)
accuracy(train.lm6.pred, valid.ts)
plot(train.lm6.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

#Building linear model with multiplicative seasonality and trend
train.lm4 <- tslm(train.ts ~ season * trend)
train.lm4.pred <- forecast(train.lm4, h = nValid)
accuracy(train.lm4.pred, valid.ts)
plot(train.lm4.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

# Building linear model with multiplicative seasonality and exponential trend
train.lm5 <- tslm(train.ts ~ season * (trend + I(trend^2)))
train.lm5.pred <- forecast(train.lm5, h = nValid)
accuracy(train.lm5.pred, valid.ts)
plot(train.lm5.pred)  
lines(valid.ts, col = "grey20", lty = 3) 

# Building Exponential models 
expo_model <- stlf(train.ts)
expo_mode_forecast <- forecast(expo_model, h = nValid)
accuracy(expo_mode_forecast, valid.ts)
plot(expo_mode_forecast)  
lines(valid.ts, col = "grey20", lty = 3) 


# So far we can see that the best model is the linear model with seasonality 
# We will test to see if there is auto-correlation
season_data <- as.data.frame(train.lm1.pred)
valid_data <- as.data.frame(valid.ts)

error <- season_data[,1] - valid_data
Acf(error)

# We can see that there is auto-correlation from lags 1-5, meaning the values 
# observed for a given month is being impacted by the previous 5 months

# Using ARIMA model to deal with the Auto-correlation
arima_model <- auto.arima(train.ts)
arima_model_forecast <- forecast(arima_model, h = nValid)
accuracy(arima_model_forecast, valid.ts)
plot(arima_model_forecast)  
lines(valid.ts, col = "grey20", lty = 3) 

# Multiple plots of ARIMA since it seems like a pretty good model
# Plot of only the forecasted part (with CI) with the actual data in the background 
plot(arima_model_forecast, xlim = c(2021.5,2022.5), ylab = "PM 2.5", xlab = "Years")  
lines(valid.ts, col = "grey20", lty = 3) 

# Plot all the forecasted part (with CI) with the actual data in the background 
plot(arima_model_forecast, ylab = "PM 2.5", xlab = "Years")  
lines(valid.ts, col = "grey20", lty = 3) 
lines(arima_model_forecast$fitted, col = "blue", lwd = 2)

# Plot the residuals of the forecast
plot(arima_model_forecast$residuals, ylab = "Forecast Errors", xlab = "Years")

#-------------------------------------------------------------------------------

# The best model is still the Seasonal Naive Model, 
# despite it not accounting for Auto-correlation. This is since the MAE 
# values are closer between the validation and training sets. Therefore,
# we are going to re-build it but using all data we have (valid.ts & train.ts)
train.snaive.final <- snaive(weekly, h = nValid)
train.pred.snaive.final <- forecast(train.snaive.final, h = nValid)
plot(train.pred.snaive.final)  


# If we decide that dealing with auto-correlation is too important to leave to 
# the side, then we can use ARIMA on full model
arima_model_full <- auto.arima(weekly)
arima_model_forecast_full <- forecast(arima_model_full, h = nValid)

# Plotting full forecast and previous data
plot(arima_model_forecast_full)  

# Plotting only the forecasted data
plot(arima_model_forecast_full, xlim = c(2022.5,2023.4), ylab = "PM 2.5", xlab = "Years", main = "ARIMA Forecast") 


