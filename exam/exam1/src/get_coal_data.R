# script to write coal data from boot package
df  <- boot::coal
# write to txt file
write.table(
    df, file = "exam/data/coal_data.txt",
    row.names = FALSE, col.names = TRUE,
    sep = "\t")
