
source("mpcs.R")

# data
csv_names <- c("delusional_disorder.csv",
               "schizoaffective_disorder.csv",
               "schizophrenia.csv",
               "schizophreniform_disorder.csv",
               "speech_sound_disorder.csv")
disorder_list <- vector(mode="list", length=length(csv_names))

# other data
# "major_depressive_disorder_sample.csv",
# "persistent_depressive_disorder_sample.csv",
# "panic_disorder_sample.csv",
# "generalized_anxiety_disorder_sample.csv",

# reading the CSVs from the subfolder "disorders"
for(i in 1:length(csv_names)){
  disorder_list[[i]] <- read.csv(paste0("./disorders/", csv_names[i]), header=T)
}

# transforming and combining them
df_all_disorder <- plyr::join_all(disorder_list,type="full",match="first")
df_all_disorder[is.na(df_all_disorder)] <- 0
df_all_disorder_names <- unique(df_all_disorder$name)
sim_matrix <- diag(1, length((df_all_disorder_names)))

# test
# k <- 0, progress iterator
for(i in 1:(length(df_all_disorder_names)-1)){
  for(j in (i+1):length(df_all_disorder_names)){
    selected_disorder1 <- df_all_disorder_names[i]
    sd1 <- df_all_disorder[df_all_disorder$name==selected_disorder1,]
    sd1 <- sd1[,c(T,colSums(sd1[!names(df_all_disorder)=="name"])>0)]
    
    selected_disorder2 <- df_all_disorder_names[j]
    sd2 <- df_all_disorder[df_all_disorder$name==selected_disorder2,]
    sd2 <- sd2[,c(T,colSums(sd2[!names(df_all_disorder)=="name"])>0)]
    
    merged <- merge(sd1, sd2, all=T)
    merged[is.na(merged)] <- 0
    
    m_sd1 <- merged[merged$name==selected_disorder1,][,!names(merged)=="name"]
    m_sd2 <- merged[merged$name==selected_disorder2,][,!names(merged)=="name"]
    m_sd1 <- as.matrix(m_sd1)
    m_sd2 <- as.matrix(m_sd2)
    colnames(m_sd1) <- c()
    colnames(m_sd2) <- c()
    
    sim_matrix[j,i] <- round(MPCS(m_sd1, m_sd2, agg="mean"), 3)
    #k <- k + 1
  }
  # loop progress
  # print(paste0( round(k/sum(1:(length(df_all_disorder_names)-1)),2)*100 ,"% done"))
}

colnames(sim_matrix) <- df_all_disorder_names
rownames(sim_matrix) <- df_all_disorder_names

sim_matrix[upper.tri(sim_matrix, diag=F)] <- NA
show(sim_matrix)
