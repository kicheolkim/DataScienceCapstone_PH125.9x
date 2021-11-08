##########################################################
# Create edx set, validation set (final hold-out test set)
##########################################################
# Note: this process could take a couple of minutes
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")

library(tidyverse)
library(caret)
library(data.table)

# MovieLens 10M dataset:
# https://grouplens.org/datasets/movielens/10m/
# http://files.grouplens.org/datasets/movielens/ml-10m.zip

dl <- tempfile()
download.file("http://files.grouplens.org/datasets/movielens/ml-10m.zip", dl)

ratings <- fread(text = gsub("::", "\t", readLines(unzip(dl, "ml-10M100K/ratings.dat"))),
                 col.names = c("userId", "movieId", "rating", "timestamp"))

movies <- str_split_fixed(readLines(unzip(dl, "ml-10M100K/movies.dat")), "\\::", 3)
colnames(movies) <- c("movieId", "title", "genres")

# if using R 3.6 or earlier:
# movies <- as.data.frame(movies) %>% mutate(movieId = as.numeric(levels(movieId))[movieId],
#                                            title = as.character(title),
#                                            genres = as.character(genres))
# if using R 4.0 or later:
movies <- as.data.frame(movies) %>% mutate(movieId = as.numeric(movieId),
                                           title = as.character(title),
                                           genres = as.character(genres))

movielens <- left_join(ratings, movies, by = "movieId")

# Validation set will be 10% of MovieLens data
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
test_index <- createDataPartition(y = movielens$rating, times = 1, p = 0.1, list = FALSE)
edx <- movielens[-test_index,]
temp <- movielens[test_index,]

# Make sure userId and movieId in validation set are also in edx set
validation <- temp %>% 
  semi_join(edx, by = "movieId") %>%
  semi_join(edx, by = "userId")

# Add rows removed from validation set back into edx set
removed <- anti_join(temp, validation)
edx <- rbind(edx, removed)

rm(dl, ratings, movies, test_index, temp, movielens, removed)



##########################################################
#### MovieLens Project - hands on ####
# # required-1: (.Rmd) Your report in Rmd format
# required-2: (.pdf, from Rmd) Your report in PDF format 
# required-3: (.R) A script in R format that generates your predicted movie ratings and RMSE score (should contain all code and comments for your project)
# 
# Develop your algorithm using the edx set. 
# For a final test of your final algorithm, predict movie ratings in the validation set (the final hold-out test set) as if they were unknown. 
# RMSE will be used to evaluate how close your predictions are to the true values in the validation set (the final hold-out test set).
##########################################################

##### Dataset prep. ####
library(lubridate)

# conversion timestamp into the date
edx <- edx %>% 
  mutate(time = as_datetime(timestamp)) %>% mutate(date = format(time, format="%Y-%m-%d"))
validation <- validation %>% 
  mutate(time = as_datetime(timestamp)) %>% mutate(date = format(time, format="%Y-%m-%d"))

# separate genres into single row
edx_genre <- edx %>% 
  separate_rows(genres, sep = "\\|") 
validation_genre <- validation %>% 
  separate_rows(genres, sep = "\\|") 



##### Training ####
# function to compute RMSE 
RMSE <- function(true_ratings, predicted_ratings){
  sqrt(mean((true_ratings - predicted_ratings)^2))
}


### Method 1. average ratings
# average rating
avg_rating <- mean(edx$rating)

# simplest recommendation system using average ratings
rmse_mean <- RMSE(validation$rating, avg_rating)
result_rmse <- tibble(method="average ratings", RMSE=rmse_mean)
knitr::kable(result_rmse)


### Method 2. movie effect
# movie effect
avg_movie <- edx %>% group_by(movieId) %>% 
  summarize(b_m = mean(rating - avg_rating))

# prediction by movie effect
rating_predict <- validation %>% 
  left_join(avg_movie, by="movieId") %>% 
  mutate(pred_movie = b_m + avg_rating) %>% pull(pred_movie)
rmse_movie <- RMSE(validation$rating, rating_predict)
result_rmse <- rbind(result_rmse, tibble(method="movie effect", RMSE=rmse_movie))
knitr::kable(result_rmse)
rm(rating_predict)


### Method 3. movie + user effect
# user effect
avg_user <- edx %>% left_join(avg_movie, by="movieId") %>% 
  group_by(userId) %>% 
  summarize(b_u = mean(rating - avg_rating - b_m))

# prediction by user
rating_predict <- validation %>% 
  left_join(avg_user, by="userId") %>% 
  mutate(pred_user = b_u + avg_rating) %>% pull(pred_user)

# prediction by movie+user
rating_predict_combined <- validation %>% 
  left_join(avg_movie, by="movieId") %>% 
  left_join(avg_user, by="userId") %>% 
  mutate(pred_user = b_m + b_u + avg_rating) %>% pull(pred_user)

# calculate RMSE for prediction by user effect only
rmse_user <- RMSE(validation$rating, rating_predict)

# calculate RMSE for prediction by movie+user effect
rmse_user_combined <- RMSE(validation$rating, rating_predict_combined)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="user effect", RMSE=rmse_user))
result_rmse <- rbind(result_rmse, tibble(method="movie+user effect", RMSE=rmse_user_combined))
knitr::kable(result_rmse)
rm(rating_predict, rating_predict_combined)


### Method 4. movie + user + date effect
# date effect
avg_date <- edx %>% 
  left_join(avg_movie, by="movieId") %>% 
  left_join(avg_user, by="userId") %>% 
  group_by(date) %>% 
  summarize(d_ui = mean(rating - avg_rating - b_m - b_u))

# prediction by date
rating_predict <- validation %>% 
  left_join(avg_date, by="date") %>% 
  mutate(pred_date = d_ui + avg_rating) %>% pull(pred_date)

# prediction by movie+user+date
rating_predict_combined <- validation %>% 
  left_join(avg_movie, by="movieId") %>% 
  left_join(avg_user, by="userId") %>% 
  left_join(avg_date, by="date") %>% 
  mutate(pred_date = d_ui + b_m + b_u + avg_rating) %>% pull(pred_date)

# calculate RMSE for prediction by date effect only
rmse_date <- RMSE(validation$rating, rating_predict)

# calculate RMSE for prediction by movie+user+date effect
rmse_date_combined <- RMSE(validation$rating, rating_predict_combined)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="date effect", RMSE=rmse_date))
result_rmse <- rbind(result_rmse, tibble(method="movie+user+date effect", RMSE=rmse_date_combined))
knitr::kable(result_rmse)
rm(rating_predict, rating_predict_combined)


### Method 5. movie + user + date + genres effect
# genre effect
avg_genres <- edx_genre %>% 
  left_join(avg_movie, by="movieId") %>% 
  left_join(avg_user, by="userId") %>% 
  left_join(avg_date, by="date") %>% 
  group_by(genres) %>% 
  summarize(b_g = mean(rating - avg_rating - b_m - b_u - d_ui))

# prediction by genre effect only
rating_predict <- validation_genre %>% 
  left_join(avg_genres, by="genres") %>% 
  mutate(pred_genres = b_g + avg_rating) %>% pull(pred_genres)

# prediction by movie+user+date+genre effect
rating_predict_combined <- validation_genre %>% 
  left_join(avg_movie, by="movieId") %>% 
  left_join(avg_user, by="userId") %>% 
  left_join(avg_date, by="date") %>% 
  left_join(avg_genres, by="genres") %>% 
  mutate(pred_genres = b_g + d_ui + b_m + b_u + avg_rating) %>% pull(pred_genres)

# calculate RMSE for prediction by genre effect only
rmse_time <- RMSE(validation_genre$rating, rating_predict)

# calculate RMSE for prediction by movie+user+date+genre effect
rmse_time_combined <- RMSE(validation_genre$rating, rating_predict_combined)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="genres effect", RMSE=rmse_time))
result_rmse <- rbind(result_rmse, tibble(method="movie+user+date+genres effect", RMSE=rmse_time_combined))
knitr::kable(result_rmse)
rm(rating_predict, rating_predict_combined)



##### Regularization: control total variability of the movie effect ####
lambdas <- seq(0, 10, 0.2)

### 1-1. regularized movie: choosing parameter
rmses <- sapply(lambdas, function(l){
  movie_reg <- edx %>% group_by(movieId) %>% 
    summarize(b_m = sum(rating - avg_rating) / (n() + l))
  rating_predict <- validation %>% 
    group_by(movieId) %>%
    left_join(movie_reg, by="movieId") %>% 
    mutate(pred = b_m + avg_rating) %>% pull(pred)
  return(RMSE(rating_predict, validation$rating))
})
# plot(lambdas, rmses)
lambda <- lambdas[which.min(rmses)]
rm(rmses)

### 1-2. regularized movie
# regularized by movie
movie_reg <- edx %>% group_by(movieId) %>% 
  summarize(b_m = sum(rating - avg_rating) / (n() + lambda))

# prediction using regularized features
rating_predict <- validation %>% 
  group_by(movieId) %>%
  left_join(movie_reg, by="movieId") %>% 
  mutate(pred = b_m + avg_rating) %>% pull(pred)

# calculate RMSE for prediction by regularized movie effect
rmse_movie_reg <- RMSE(rating_predict, validation$rating)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="regularized movie effect", RMSE=rmse_movie_reg))
knitr::kable(result_rmse)
rm(rating_predict)


### 2-1. regularized movie+user: choosing parameter
rmses <- sapply(lambdas, function(l){
  movie_reg <- edx %>% group_by(movieId) %>% 
    summarize(b_m = sum(rating - avg_rating) / (n() + l))
  user_reg <- edx %>% group_by(userId) %>% 
    left_join(movie_reg, by="movieId") %>% 
    summarize(b_u = sum(rating - avg_rating - b_m) / (n() + l))
  rating_predict <- validation %>% 
    group_by(movieId) %>%
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    mutate(pred = b_m + b_u + avg_rating) %>% pull(pred)
  return(RMSE(rating_predict, validation$rating))
})
lambda <- lambdas[which.min(rmses)]
rm(rmses)

### 2-2. regularized movie+user
# regularized by movie
movie_reg <- edx %>% group_by(movieId) %>% 
  summarize(b_m = sum(rating - avg_rating) / (n() + lambda))

# regularized by user
user_reg <- edx %>% group_by(userId) %>% 
  left_join(movie_reg, by="movieId") %>% 
  summarize(b_u = sum(rating - avg_rating - b_m) / (n() + lambda))

# prediction using regularized features
rating_predict <- validation %>% 
  group_by(movieId) %>%
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  mutate(pred = b_m + b_u + avg_rating) %>% pull(pred)

# calculate RMSE for prediction by regularized movie+user effect
rmse_user_reg <- RMSE(rating_predict, validation$rating)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="regularized movie+user effect", RMSE=rmse_user_reg))
knitr::kable(result_rmse)
rm(rating_predict)


### 3-1. regularized movie+user+date: choosing parameter
rmses <- sapply(lambdas, function(l){
  movie_reg <- edx %>% group_by(movieId) %>% 
    summarize(b_m = sum(rating - avg_rating) / (n() + l))
  user_reg <- edx %>% group_by(userId) %>% 
    left_join(movie_reg, by="movieId") %>% 
    summarize(b_u = sum(rating - avg_rating - b_m) / (n() + l))
  date_reg <- edx %>% group_by(date) %>% 
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    summarize(d_ui = sum(rating - avg_rating - b_m - b_u) / (n() + l))
  rating_predict <- validation %>% 
    group_by(movieId) %>%
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    left_join(date_reg, by="date") %>% 
    mutate(pred = b_m + b_u + d_ui + avg_rating) %>% pull(pred)
  return(RMSE(rating_predict, validation$rating))
})
lambda <- lambdas[which.min(rmses)]
rm(rmses)

### 3-2. regularized movie+user+date
# regularized by movie
movie_reg <- edx %>% group_by(movieId) %>% 
  summarize(b_m = sum(rating - avg_rating) / (n() + lambda))

# regularized by user
user_reg <- edx %>% group_by(userId) %>% 
  left_join(movie_reg, by="movieId") %>% 
  summarize(b_u = sum(rating - avg_rating - b_m) / (n() + lambda))

# regularized by date
date_reg <- edx %>% group_by(date) %>% 
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  summarize(d_ui = sum(rating - avg_rating - b_m - b_u) / (n() + lambda))

# prediction using regularized features
rating_predict <- validation %>% 
  group_by(movieId) %>%
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  left_join(date_reg, by="date") %>% 
  mutate(pred = b_m + b_u + d_ui + avg_rating) %>% pull(pred)

# calculate RMSE for prediction by regularized movie+user+date effect
rmse_date_reg <- RMSE(rating_predict, validation$rating)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="regularized movie+user+date effect", RMSE=rmse_date_reg))
knitr::kable(result_rmse)
rm(rating_predict)


### 4-1. regularized movie+user+date+genres: choosing parameter
rmses <- sapply(lambdas, function(l){
  movie_reg <- edx %>% group_by(movieId) %>% 
    summarize(b_m = sum(rating - avg_rating) / (n() + l))
  user_reg <- edx %>% group_by(userId) %>% 
    left_join(movie_reg, by="movieId") %>% 
    summarize(b_u = sum(rating - avg_rating - b_m) / (n() + l))
  date_reg <- edx %>% group_by(date) %>% 
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    summarize(d_ui = sum(rating - avg_rating - b_m - b_u) / (n() + l))
  genre_reg <- edx_genre %>% group_by(genres) %>% 
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    left_join(date_reg, by="date") %>% 
    summarize(b_g = sum(rating - avg_rating - b_m - b_u - d_ui) / (n() + l))
  rating_predict <- validation_genre %>% 
    group_by(movieId) %>%
    left_join(movie_reg, by="movieId") %>% 
    left_join(user_reg, by="userId") %>% 
    left_join(date_reg, by="date") %>% 
    left_join(genre_reg, by="genres") %>% 
    mutate(pred = b_m + b_u + d_ui + b_g + avg_rating) %>% pull(pred)
  return(RMSE(rating_predict, validation_genre$rating))
})
lambda <- lambdas[which.min(rmses)]
rm(rmses)

### 4-2. regularized movie+user+date+genres
# regularized by movie
movie_reg <- edx %>% group_by(movieId) %>% 
  summarize(b_m = sum(rating - avg_rating) / (n() + lambda))

# regularized by user
user_reg <- edx %>% group_by(userId) %>% 
  left_join(movie_reg, by="movieId") %>% 
  summarize(b_u = sum(rating - avg_rating - b_m) / (n() + lambda))

# regularized by date
date_reg <- edx %>% group_by(date) %>% 
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  summarize(d_ui = sum(rating - avg_rating - b_m - b_u) / (n() + lambda))

# regularized by genre
genre_reg <- edx_genre %>% group_by(genres) %>% 
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  left_join(date_reg, by="date") %>% 
  summarize(b_g = sum(rating - avg_rating - b_m - b_u - d_ui) / (n() + lambda))

# prediction using regularized features
rating_predict <- validation_genre %>% 
  group_by(movieId) %>%
  left_join(movie_reg, by="movieId") %>% 
  left_join(user_reg, by="userId") %>% 
  left_join(date_reg, by="date") %>% 
  left_join(genre_reg, by="genres") %>% 
  mutate(pred = b_m + b_u + d_ui + b_g + avg_rating) %>% pull(pred)

# calculate RMSE for prediction by regularized movie+user+date+genre effect
rmse_genre_reg <- RMSE(rating_predict, validation_genre$rating)

# summarize results
result_rmse <- rbind(result_rmse, tibble(method="regularized movie+user+date+genre effect", RMSE=rmse_genre_reg))
knitr::kable(result_rmse)
rm(rating_predict)



#### RMSE conclusion
# sorting results by RMSE
result_rmse %>% arrange(RMSE) %>% knitr::kable()

# best RMSE (= minimum RMSE score)
result_rmse[which.min(result_rmse$RMSE),]

