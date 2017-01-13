# Home for general utility functions

mkdir <- function(baseDir_v, 
                  newDir_v){
  # Create new file for writing (if doesn't exist)
  # baseDir_v: some sort of file path to where new files should be written
  # newDir_v:  character vector of name of new directory to write files to
  tempDir_v <- file.path(baseDir_v, newDir_v)
  if(dir.exists(tempDir_v)){
    return(tempDir_v)
  } else {
    dir.create(tempDir_v)
    return(tempDir_v)
  } # fi
} # mkdir