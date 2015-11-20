learnIcelandic <- function(csvDirPath){
  require(tm)
  require(stylo)
  
  # Get stuff from althingi_tagged, assuming most of it is correct:
  filePaths <- DirSource(csvDirPath)$filelist
  
  words <- c()
  lemmas <- c()
  for(fileName in filePaths){
    print(paste("Reading file",fileName,"..."))
    rawData <- read.csv(fileName, fileEncoding="UTF-8")
    print("done!")
    
    print("Gathering words...")
    trim1 <- rawData$Word[!grepl("[^[:alpha:]]",rawData$Word)]
    trim2 <- trim1[trim1!=""]
    words <- c(words,tolower(as.character(trim2)))
    print("done!")
    
    
    print("Gathering lemmas...")
    trim1 <- rawData$Lemma[!grepl("[^[:alpha:]]",rawData$Lemma)]
    trim2 <- trim1[trim1!=""]
    lemmas <- c(lemmas,tolower(as.character(trim2)))
    print("done!")
  }
  
  dict <<- table(words) # Get a unique list of words and their frequency, alphabetized
  # Lemma bigrams:
  bigrams <- make.ngrams(lemmas,ngram.size = 2)
  # ALL bigrams: (all words) -- doesn't work?
  word.bigrams <- make.ngrams(words,ngram.size = 2)
  word.bigramDict <<- table(word.bigrams)
  
  # Lemma bigram dictionary, with frequency
  bigramDict <<- table(bigrams)
  # Lemma dictionary, with frequency
  lemmasDict <<- table(lemmas)
}

#
find = function(thing, dict) {
  return (length(dict[grepl(paste("^",thing,"$", sep = ""), names(dict))]))
}

# Function that calculates the probability of a lemma
P_lemma = function(l) {
  freq <- as.numeric(lemmasDict[l])+1 #Frequency of lemma 'l' in the dictionary.
  if(is.na(freq)) freq <- 1 # set as 1 if it doesn't appear anywhere
  
  # Calculate the probability of the word:
  return (freq/as.numeric(max(lemmasDict)+1))
}

# Function that calculates the probability of a word
P_word = function(w) {
  freq <- as.numeric(dict[w])+1 #Frequency of lemma 'l' in the dictionary.
  if(is.na(freq)) freq <- 1 # set as 1 if it doesn't appear anywhere
  
  # Calculate the probability of the word:
  return (freq/as.numeric(max(dict)+1))
}

# Function that calculates the probability that a 
# B comes after A in a sentence
P_bigram_lemma_AX <- function(A, B) {
  bigram <- paste(A,B)
  freq <<- as.numeric(bigramDict[bigram])+1
  if(is.na(freq)) freq <- 1
  return(max(min(freq/10,0.9),0.1))
}

# Function that calculates the probability that a 
# B comes after A in a sentence
P_bigram_AX <- function(A, B) {
  bigram <- paste(A,B)
  freq <<- as.numeric(wordBigramDict[bigram])+1
  if(is.na(freq)) freq <- 1
  freqAX <<- as.numeric(sum(wordBigramDict[grepl(paste("^",A,sep=""),names(wordBigramDict))]))+1
  if(is.na(freqAX)) freqAX <- 2
  return(freq/freqAX)
}

# Function that calculates the probability that a 
# B comes after A in a sentence
P_bigram_XC <- function(B, C) {
  bigram <- paste(B,C)
  freq <- as.numeric(wordBigramDict[bigram])+1
  if(is.na(freq)) freq <- 1
  freqXC <- as.numeric(sum(wordBigramDict[grepl(paste(C,"$",sep=""),names(wordBigramDict))]))+1
  if(is.na(freqXC)) freqXC <- 2
  return(freq/freqXC)
}


# Function that calculates the normalized probability that B comes after A
P_trigram_norm <- function(A, B, C, verbose=FALSE){
  print("Calculating...")
  pAB <- P_bigram_AX(A,B)
  if(verbose) print(paste("Probability of B after A:",pAB))
  pBC <- P_bigram_XC(B,C)
  if(verbose) print(paste("Probability of B before C:",pBC))
  pCorr <- pAB*pBC
  if(verbose) print(paste("Probability the trigram is correct:",pCorr))
  pIncorr <- (1-pAB)*(1-pBC)
  if(verbose) print(paste("Probability the trigram is incorrect:",pIncorr))
  pNorm <- pCorr/(pCorr+pIncorr)
  return(pNorm)
}

# Function that calculates the normalized probability that B comes after A
P_bigram_norm <- function(A, B, verbose=FALSE){
  print("Calculating...")
  pB <- max(min(lemmasDict[B]/10,0.9),0.1)
  if(is.na(pB)) pB <- 0.1
  if(verbose) print(paste("Probability of B :",pB))
  pAB <- P_bigram_lemma_AX(A,B)
  if(verbose) print(paste("Probability of B after A:",pAB))
  pCorr <- pB * pAB
  pIncorr <- (1-pB)*(1-pAB)
  pNorm <- pCorr/(pCorr+pIncorr)
  return(pNorm)
}

# Function that takes in a string "A B C"
# and evaluates the probability that "B" makes sense
# in its' context.
naiveBayesTrigram <- function(trigram){
  require(stylo)
  # Make trigram array:
  words <- txt.to.words(trigram)
  if(length(words) != 3){
    print("Not a valid trigram!")
    return()
  }
  
  # "A B"
  frontBigram <- paste(words[1],words[2])
  # "B C"
  backBigram <- paste(words[2],words[3])
  
  # Calculate the various frequencies
  # ------------- ATH !!!!!!!!!!!!!! ---------------------------
  # bigramDict[bigram] tekur 12 missisippi á meðan
  # length(bigramDict[grepl(bigram, names(bigramDict))])
  # tekur 2 missisippi...stórfelld áhrif á hraðvirkni forritsins
  # ------------------------------------------------------------
  fbFreq <- bigramDict[frontBigram]  # "A B"
  bbFreq <- bigramDict[backBigram] # "B C"
  w1Freq <- lemmasDict[words[1]] # "A"
  w2Freq <- lemmasDict[words[2]] # "B"
  w3Freq <- lemmasDict[words[3]] # "C"
  
  # Initialize in case we get NA values:
  if(is.na(fbFreq)) fbFreq <- 1
  if(is.na(bbFreq)) bbFreq <- 1
  if(is.na(w1Freq)) w1Freq <- 1
  if(is.na(w2Freq)) w2Freq <- 1
  if(is.na(w3Freq)) w3Freq <- 1
  
  # Calculate the frequencies like so:
  # A bigram is considered to have 100% chance of being
  # correct if it appears in our bigram dictionary more 
  # than a 100 times. Otherwise it's frequency is divided by 100.
  fbP <- min(fbFreq,100)/100
  bbP <- min(bbFreq,100)/100
  w1P <- 1
  w2P <- min(w2Freq,100)/100
  w3P <- 1
  
  # Naive Bayes allows us to multiply 
  # the odds as if they were independent:
  # P(B | A) = P(A B)*P(A)*P(B)
  fbC <- fbP*w1P*w2P
  bbC <- bbP*w2P*w3P
  
  totalC <- fbC*bbC
  
  return(as.numeric(totalC))
}

# Use Lehvenstein distance to find suggestions
# for words that are deemed incorrect
suggestWord <- function(word){
  require(stringdist)
  dist <- stringdist(names(dict),word)
  possibleWords <- dict[ 0 < dist & dist < 3]
  return(sort(possibleWords,decreasing = TRUE))
}

# Select the suggestion that is likely to be correct (i.e. is common):
selectSuggestion <- function(pw, word){
  require(stringdist)
  adjust <- function(dist) return(1.7-0.8*dist)
  for(suggestion in names(pw)){
    pw[suggestion] <- pw[suggestion]*adjust(stringdist(suggestion,word))
  }
  return(names(sort(pw,decreasing = TRUE))[1])
}


# Takes a word/tag/lemma csv, reads the words and makes
# a decision for each word whether it's correct in it's context
# If a word is deemed incorrect, it put's a suggestion in it's place.
correctFile <- function(csvFile,destFile){
  require(stylo)
  rawData <- read.csv(csvFile)
  file_words <- as.character(rawData$Word)
  file_tags <- as.character(rawData$Tag)
  file_lemmas_raw <- as.character(rawData$Lemma)
  file_lemmas <- tolower(file_lemmas_raw)
  
  ##First we try to correct non-word errors in the file:
  print("Non-word error pass...")
  possibleWord <- c()
  for(i in 1:length(file_words)){
    word <- tolower(file_words[i])
    freq <- dict[word]
    correction <- word
    if(is.na(freq)){
      suggestions <- suggestWord(word)
      correction <- selectSuggestion(suggestions, word)
      print(paste(word,"is not a word, replacing with best guess:",correction))
    }
    possibleWord[i] <- correction
  }
  
  ##Then we do a context based correction
  print("Context error pass...")
  correctWord <- c(possibleWord[1])
  for(i in 2:(length(possibleWord)-1)){
    print("Working...")
    print(paste("calculating probability for:",paste(file_lemmas[i-1],file_lemmas[i])))
    probCorrect <- P_bigram_norm(file_lemmas[i-1],file_lemmas[i],verbose = TRUE)
    print(paste("result:",probCorrect))
    verdict <- probCorrect>0.1
#    print(paste("This word is correct:",verdict))
    correction <- probCorrect #possibleWord[i]
#     if(!verdict){
#       suggestions <- suggestWord(possibleWord[i])
#       correction <- selectSuggestion(suggestions, possibleWord[i])
#     } 
#     print(correction)
    correctWord <- c(correctWord,correction)
  }
  correctWord <- c(correctWord,possibleWord[length(possibleWord)])
  check <<- correctWord
  correctedData <- data.frame(Word=file_words,Tag=file_tags,Lemma=file_lemmas_raw,CorrectWord=correctWord)
  write.csv(correctedData, file=destFile)
}
