## Order function
# changes orders of columns in data frame, knows before, after, first, last.


#usage
#Here's some sample data representing the names of your dataset:
#
#x <- paste0("v", 1:18)
#Imagine now that we wanted "v17" and "v18" before "v3", "v6" and "v16" at the end, and "v5" at the beginning:
#
#moveme(x, "v17, v18 before v3; v6, v16 last; v5 first")
##  [1] "v5"  "v1"  "v2"  "v17" "v18" "v3"  "v4"  "v7"  "v8"  "v9"  "v10" "v11" "v12"
## [14] "v13" "v14" "v15" "v6"  "v16"

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}     