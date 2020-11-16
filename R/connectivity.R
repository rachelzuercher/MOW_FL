# calcuating larval connectivity for FL reef tract
# Mapping Ocean Wealth FL

library(R.matlab)

locations <- read.table("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/polyKeys10k.xyz")
colnames(locations)[1] <- "lon"
colnames(locations)[2] <- "lat"
colnames(locations)[3] <- "id"

# 2005 #
mat2005_4 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month4.mat")
mat2005_4 <- as.matrix(mat2005_4$settle) # extract just the data matrix from the Matlab list
mat2005_4 <- rbind(mat2005_4,NA,NA) # add two rows for summing the recruitment data
colnames(mat2005_4) <- as.character(1:50)
rownames(mat2005_4) <- as.character(1:52)
mat2005_5 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month5.mat")
mat2005_5 <- as.matrix(mat2005_5$settle)
mat2005_5 <- rbind(mat2005_5,NA,NA)
colnames(mat2005_5) <- as.character(1:50)
rownames(mat2005_5) <- as.character(1:52)
mat2005_6 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month6.mat")
mat2005_6 <- as.matrix(mat2005_6$settle)
mat2005_6 <- rbind(mat2005_6,NA,NA)
colnames(mat2005_6) <- as.character(1:50)
rownames(mat2005_6) <- as.character(1:52)
mat2005_7 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month7.mat")
mat2005_7 <- as.matrix(mat2005_7$settle)
mat2005_7 <- rbind(mat2005_7,NA,NA)
colnames(mat2005_7) <- as.character(1:50)
rownames(mat2005_7) <- as.character(1:52)
mat2005_8 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month8.mat")
mat2005_8 <- as.matrix(mat2005_8$settle)
mat2005_8 <- rbind(mat2005_8,NA,NA)
colnames(mat2005_8) <- as.character(1:50)
rownames(mat2005_8) <- as.character(1:52)
mat2005_9 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month9.mat")
mat2005_9 <- as.matrix(mat2005_9$settle)
mat2005_9 <- rbind(mat2005_9,NA,NA)
colnames(mat2005_9) <- as.character(1:50)
rownames(mat2005_9) <- as.character(1:52)
mat2005_10 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month10.mat")
mat2005_10 <- as.matrix(mat2005_10$settle)
mat2005_10 <- rbind(mat2005_10,NA,NA)
colnames(mat2005_10) <- as.character(1:50)
rownames(mat2005_10) <- as.character(1:52)
mat2005_11 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2005_month11.mat")
mat2005_11 <- as.matrix(mat2005_11$settle)
mat2005_11 <- rbind(mat2005_11,NA,NA)
colnames(mat2005_11) <- as.character(1:50)
rownames(mat2005_11) <- as.character(1:52)

array2005 <- array(c(mat2005_4, mat2005_5, mat2005_6, mat2005_7, mat2005_8, mat2005_9, mat2005_10, mat2005_11), dim=c(52,50,8))
rm(mat2005_4, mat2005_5, mat2005_6, mat2005_7, mat2005_8, mat2005_9, mat2005_10, mat2005_11)

# column totals (total recruitment in Row 51)
for (i in 1:8) {
  array2005[51,,i] <- colSums(array2005[1:50,,i])
}

# row 52 is col sums minus self-recruitment
for (i in 1:8) {
  for (j in 1:50) {
    array2005[52,j,i] <- array2005[51,j,i]-array2005[j,j,i]
  }
}

# sum all months of a given year for each location 1-50 (i.e. add the value in row 51 for each sheet of the array for each column)
monthly_totals2005 <- data.frame(matrix(NA, nrow = 1, ncol = 50))
colnames(monthly_totals2005) <- as.character(1:50)

for (i in 1:50) {
  temp <- array2005[51,i,1] + array2005[51,i,2] + array2005[51,i,3] + array2005[51,i,4] + array2005[51,i,5] +
    array2005[51,i,6] + array2005[51,i,7] + array2005[51,i,8]
  monthly_totals2005[i] <- temp
}



# 2006 #
mat2006_4 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month4.mat")
mat2006_4 <- as.matrix(mat2006_4$settle) # extract just the data matrix from the Matlab list
mat2006_4 <- rbind(mat2006_4,NA,NA) # add two rows for summing the recruitment data
colnames(mat2006_4) <- as.character(1:50)
rownames(mat2006_4) <- as.character(1:52)
mat2006_5 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month5.mat")
mat2006_5 <- as.matrix(mat2006_5$settle)
mat2006_5 <- rbind(mat2006_5,NA,NA)
colnames(mat2006_5) <- as.character(1:50)
rownames(mat2006_5) <- as.character(1:52)
mat2006_6 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month6.mat")
mat2006_6 <- as.matrix(mat2006_6$settle)
mat2006_6 <- rbind(mat2006_6,NA,NA)
colnames(mat2006_6) <- as.character(1:50)
rownames(mat2006_6) <- as.character(1:52)
mat2006_7 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month7.mat")
mat2006_7 <- as.matrix(mat2006_7$settle)
mat2006_7 <- rbind(mat2006_7,NA,NA)
colnames(mat2006_7) <- as.character(1:50)
rownames(mat2006_7) <- as.character(1:52)
mat2006_8 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month8.mat")
mat2006_8 <- as.matrix(mat2006_8$settle)
mat2006_8 <- rbind(mat2006_8,NA,NA)
colnames(mat2006_8) <- as.character(1:50)
rownames(mat2006_8) <- as.character(1:52)
mat2006_9 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month9.mat")
mat2006_9 <- as.matrix(mat2006_9$settle)
mat2006_9 <- rbind(mat2006_9,NA,NA)
colnames(mat2006_9) <- as.character(1:50)
rownames(mat2006_9) <- as.character(1:52)
mat2006_10 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month10.mat")
mat2006_10 <- as.matrix(mat2006_10$settle)
mat2006_10 <- rbind(mat2006_10,NA,NA)
colnames(mat2006_10) <- as.character(1:50)
rownames(mat2006_10) <- as.character(1:52)
mat2006_11 <- readMat("/Users/rachelzuercher/Desktop/FIU/TNC_Florida/Data/Biophysical/connectivity/matrices/new_matrix_raw_year2006_month11.mat")
mat2006_11 <- as.matrix(mat2006_11$settle)
mat2006_11 <- rbind(mat2006_11,NA,NA)
colnames(mat2006_11) <- as.character(1:50)
rownames(mat2006_11) <- as.character(1:52)

array2006 <- array(c(mat2006_4, mat2006_5, mat2006_6, mat2006_7, mat2006_8, mat2006_9, mat2006_10, mat2006_11), dim=c(52,50,8))
rm(mat2006_4, mat2006_5, mat2006_6, mat2006_7, mat2006_8, mat2006_9, mat2006_10, mat2006_11)

# column totals (total recruitment in Row 51)
for (i in 1:8) {
  array2006[51,,i] <- colSums(array2006[1:50,,i])
}

# row 52 is col sums minus self-recruitment
for (i in 1:8) {
  for (j in 1:50) {
    array2006[52,j,i] <- array2006[51,j,i]-array2006[j,j,i]
  }
}

# sum all months of a given year for each location 1-50 (i.e. add the value in row 51 for each sheet of the array for each column)
monthly_totals2006 <- data.frame(matrix(NA, nrow = 1, ncol = 50))
colnames(monthly_totals2006) <- as.character(1:50)

for (i in 1:50) {
  temp <- array2006[51,i,1] + array2006[51,i,2] + array2006[51,i,3] + array2006[51,i,4] + array2006[51,i,5] +
    array2006[51,i,6] + array2006[51,i,7] + array2006[51,i,8]
  monthly_totals2006[i] <- temp
}


# average across two years #

rownames(monthly_totals2005)[1] <- "total2005"
rownames(monthly_totals2006)[1] <- "total2006"
totals <- rbind(monthly_totals2005, monthly_totals2006)
rm(monthly_totals2005, monthly_totals2006)
totals <- rbind(totals, NA)
rownames(totals)[3] <- "annual.average"
for (i in 1:50) {
  totals[3,i] <- (totals[1,i] + totals[2,i]) / 2
}

connectivity <- totals[3,]

connectivity <- rbind(connectivity, colnames(connectivity))
rownames(connectivity)[2] <- "connectivity_id"
connectivity <- as.data.frame(t(connectivity))
connectivity <- connectivity[,c(2,1)]
connectivity$annual.average <- as.numeric(as.character(connectivity$annual.average))
colnames(connectivity)[2] <- "connectivity"

# process locations data #
locations$lon <- locations$lon-360

rm(locations, totals, array2005, array2006, i, j, temp)

fwrite(connectivity, "/Users/rachelzuercher/Desktop/connectivity.csv")
