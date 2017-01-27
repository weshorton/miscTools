outputData_dt <- sapply((1:metadata_dt[,.N]), function(i){
  if (i %% 10 == 0){cat(i, " ")}
  if (i == 1){
    comboData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
  } else {
    toAddData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
    # Outer join data
    comboData_dt <- merge(x = comboData_dt, y = toAddData_dt, by = obsCol_v, all = T)
  } # fi
  return(comboData_dt)})

for (i in 1:metadata_dt[,.N]){
  if (i %% 10 == 0){cat(i)}
  if (i == 1){
    comboData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
  } else {
    toAddData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
    # Outer join data
    comboData_dt <- merge(x = comboData_dt, y = toAddData_dt, by = obsCol_v, all = T)
  } # fi
} # for