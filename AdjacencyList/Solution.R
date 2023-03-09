library(Biostrings)
library(plyr)

adjacencyList <- function(fastaFile, k){ #k = nuber of nucleotides to compare
  data <- readDNAStringSet(fastaFile)# import dataset and extract only sequences
  sequences <- matrix(data)
  sequences <- (strsplit(sequences,''))
  sequences <- plyr::ldply(sequences, rbind)
  graph <- NULL
  answer <- data.frame()
  ID <- (as.data.frame(data@ranges))$names #store sample IDs
  for (i in 1:nrow(sequences)){
    for (j in 1:nrow(sequences)){
      end_seq <-(sequences[i,
                           seq(length(sequences[i,])-
                                 sum(is.na(sequences[i,]))-(k-1),length(sequences[i,])-
                                 sum(is.na(sequences[i,])))
      ])
      
      start_seq <-sequences[j,1:k]
      if (paste(end_seq, collapse = '') == paste(start_seq, collapse = '') ){#compare nucleotides from ends of sequences to the start of others
        graph <- rbind(graph, c(i,j))
      }
    }
  }
  t = 1
  while(t <= nrow(graph)){ #this removes loops from the graph (IDs pointing to themselves)
    if (graph[t,1] == graph[t,2]){
      graph <- graph[-(t),]
    }else{
        t = t + 1
      }
    }
   
  for (i in 1:nrow(graph)){ #this substitutes the numbers by corresponding IDs
    answer[i,1] <- ID[graph[i,1]]
    answer[i,2] <- ID[graph[i,2]] 
  }
  return(answer)
}

resposta = adjacencyList('rosalind_grph.txt', 3)

write.csv(resposta, 'resposta_teste2')

