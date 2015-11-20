learnIcelandic <- function(csvDirPath,verbose=FALSE){
  require(tm)
  require(stylo)
  
  # Get stuff from althingi_tagged, assuming most of it is correct:
  filePaths <- DirSource(csvDirPath)$filelist
  
  words <- c()
  lemmas <- c()
  for(fileName in filePaths){
    if(verbose)print(paste("Reading file",fileName,"..."))
    rawData <- read.csv(fileName, fileEncoding="UTF-8")
    if(verbose)print("done!")
    
    if(verbose)print("Gathering words...")
    trim1 <- rawData$Word[!grepl("[^[:alpha:]]",rawData$Word)]
    trim2 <- trim1[trim1!=""]
    words <- c(words,tolower(as.character(trim2)))
    if(verbose)print("done!")
    
    
    if(verbose)print("Gathering lemmas...")
    trim1 <- rawData$Lemma[!grepl("[^[:alpha:]]",rawData$Lemma)]
    trim2 <- trim1[trim1!=""]
    lemmas <- c(lemmas,tolower(as.character(trim2)))
    if(verbose)print("done!")
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

# Function that calculates the probability of a lemma
P_lemma = function(lemma) {
  n <- 10
  if(nchar(lemma)==1) n <- 10000
  pL <- max(min(lemmasDict[tolower(lemma)]/n,0.9),0.1)
  if(is.na(pL)) pL <- 0.1
  
  return(pL)
}

# Function that calculates the probability that a 
# B comes after A in a sentence
P_B_after_A <- function(A, B) {
  n <- 100
  bigram <- paste(A,B)
  freq <<- as.numeric(bigramDict[bigram])
  if(is.na(freq)) return(0.1)
  return(max(min(freq/n,0.9),0.1))
}

# Function that calculates the normalized probability that B comes after A
P_bigram <- function(A, B, verbose=FALSE){
  if(verbose)print("Calculating...")
  pB <- P_lemma(B)
  if(is.na(pB)) pB <- 0.1
  if(verbose) print(paste("Probability of B :",pB))
  pAB <- P_B_after_A(A,B)
  if(verbose) print(paste("Probability of B after A:",pAB))
  pCorr <- pB * pAB
  pIncorr <- (1-pB)*(1-pAB)
  pNorm <- pCorr/(pCorr+pIncorr)
  return(pNorm)
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
  adjust <- function(dist){
    if(dist==1) return(1)
    return(0.01)
  }
  for(suggestion in names(pw)){
    pw[suggestion] <- pw[suggestion]*adjust(stringdist(suggestion,word))
  }
  return(names(sort(pw,decreasing = TRUE))[1])
}

# Give a best guess for a lemma
suggestLemma <- function(word){
  lemmaSuggestions <- wordLemmaDict[wordLemmaDict$Word==word,(2:3)]
  if(nrow(lemmaSuggestions) == 0) return(word)
  lemmaGuess <- as.character(lemmaSuggestions[lemmaSuggestions$Freq == max(lemmaSuggestions$Freq),1])
  return(lemmaGuess)
}

# Check if word is punctuation
isPunct <- function(word){
  return(grepl("[[:punct:][:digit:]]",word) | grepl("^[[:space:]]*$",word))
}

# Takes a word/tag/lemma csv, reads the words and makes
# a decision for each word whether it's correct in it's context
# If a word is deemed incorrect, it put's a suggestion in it's place.
correctFile <- function(csvFile,destFile,verbose=FALSE,extraVerbose=FALSE){
  require(stylo)
  rawData <- read.csv(csvFile)
  file_words <- as.character(rawData$Word)
  file_tags <- as.character(rawData$Tag)
  file_lemmas_raw <- as.character(rawData$Lemma)
  file_lemmas <- tolower(file_lemmas_raw)
  file_length <- length(file_words)
  
  ##First we try to correct non-word errors in the file:
  if(verbose)print("Non-word error pass...")
  possibleWord <- c()
  possibleLemma <- c()
  for(i in 1:file_length){
    progress <- floor((i*100)/file_length)
    percentage <- paste(formatC(progress, width=3, flag=" "),"%   ",sep="")
    word <- tolower(file_words[i])
    if(isPunct(word)){
      possibleLemma[i] <- file_lemmas[i]
      possibleWord[i] <- file_words[i]
    }else{
      freq <- dict[word]
      correction <- file_words[i]
      correctionLemma <- file_lemmas[i]
      if(is.na(freq)){
        suggestions <- suggestWord(word)
        correction <- selectSuggestion(suggestions, word)
        correctionLemma <- suggestLemma(correction)
        if(verbose)print(paste(percentage,word,"is not a word, replacing with best guess:",correction))
      }
      possibleLemma[i] <- correctionLemma
      possibleWord[i] <- correction
    }
  }
  
  ##Then we do a context based correction
  if(verbose)print("Context error pass...")
  startOfSentence <- TRUE
  correctWord <- c()
  firstWords <<- c()
  for(i in 1:file_length){
    progress <- floor((i*100)/file_length)
    percentage <- paste(formatC(progress, width=3, flag=" "),"%   ",sep="")
    if(startOfSentence){
      correctWord[i] <- possibleWord[i]
      firstWords <<- c(firstWords,possibleWord[i])
      startOfSentence <- FALSE
      if(verbose)print(paste(percentage,"================ NEW SENTANCE STARTING WITH",possibleWord[i],"============"))
    }else if(isPunct(possibleWord[i])){
      correctWord[i] <- possibleWord[i]
      startOfSentence <- TRUE
    }else{
      if(verbose)print(paste(percentage,"Working..."))
      if(verbose)print(paste(percentage,"calculating probability for:",paste(possibleLemma[i-1],possibleLemma[i])))
      probCorrect <- P_bigram(possibleLemma[i-1],possibleLemma[i],extraVerbose)
      if(verbose)print(paste(percentage,"result:",probCorrect))
      verdict <- floor(probCorrect*10000) > 5000
      if(verbose)print(paste(percentage,"This word is correct:",verdict))
      correction <- possibleWord[i]
      correctionLemma <- possibleLemma[i]
      if(!verdict){
        suggestions <- suggestWord(possibleWord[i])
        correction <- selectSuggestion(suggestions, possibleWord[i])
        correctionLemma <- suggestLemma(correction)
        if(is.na(correction)){
          correction <- possibleWord[i]
          correctionLemma <- possibleLemma[i]
        }
      } 
      if(verbose)print(paste(percentage,correction))
      possibleLemma[i] <- correctionLemma
      correctWord[i] <- correction
    }
  }
  check <<- correctWord
  correctedData <- data.frame(Word=file_words,Tag=file_tags,Lemma=file_lemmas_raw,CorrectWord=correctWord)
  write.csv(correctedData, file=destFile)
}

testAccuracy <- function(csvCorrectedFile, csvCorrectFile) {
  require(stylo)
  ourData <- read.csv(csvCorrectedFile)
  originalData <- read.csv(csvCorrectFile)
  ourCorrections <- ourData$CorrectWord
  orgCorrections <- originalData$CorrectWord
  
  correctCorrections <- ourCorrections[tolower(ourCorrections) == tolower(orgCorrections)]
  accuracy <- length(correctCorrections)/length(ourCorrections)
  return(accuracy)
}
