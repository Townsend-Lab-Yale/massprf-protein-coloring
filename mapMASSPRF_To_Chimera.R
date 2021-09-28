#example call
genMASSPRF_Chimera(pdbFile = "~/../Downloads/SARS-CoV-1_ORF7a_Gene.pdb",MASSPRF_Nuc_Fasta = "~/../Downloads/orf7a_gene.fasta",MASSPRF_Table = "~/../Downloads/orf7a_MASS_PRF.tsv",scaling=1,sigSetting = "average")
#genMASSPRF_Chimera(pdbFile = "~/../Desktop/S_23.pdb",MASSPRF_Nuc_Fasta = "~/../Desktop/s_gene.fasta",onlySig=F,MASSPRF_Table = "~/../Downloads/S1_MASSPRF.tsv",scaling=6,outfile = "S_cc.txt")

#pdbFile: The protein structure you are trying to map onto, in pdb format
#MASSPRF_Nuc_Fasta: the nucleotide sequence representing your MASSPRF input, in open reading frame and in fasta format
#MASSPRF_Table: The output table from MASSPRF containing the gamma values and CI
#outfile: where the resultant commands for Chimera to color will be stored
### (hint: execute the following on the Chimera command line): read <outfile>;
#scaling: the scale factor from MASSPRF--look in the output of MASSPRF to find this. 1 or a multiple of 3.
#sigSetting: what are you using to determine significance? only affects scale factor 1; otherwise sitewise significance is used
#onlySig: T or F do you want to only color significant sites
#rgb1 and rgb2: vectors of size 3, corresponding to a rgb value for coloring the selection intensity
#bins: how many equally spaced apart color categories you want to use for your data 
#ehColor: what do you want to color 'eh' residues (missing data, nonsignificant) as a vector size three of rgb


genMASSPRF_Chimera<-function(pdbFile,MASSPRF_Nuc_Fasta,MASSPRF_Table,outfile="chimeraColoring.txt",scaling=0, sigSetting="average",onlySig=T,rgb1=c(250,30,30),rgb2=c(30,30,250),bins=100,ehColor=c(128,128,128)){
  require("bio3d")
  require("Biostrings")
  require("muscle")
  if(!sigSetting%in%c("average","any","majority","strict")){stop("sigSetting must be average, any, majority, strict")}
  if(scaling==0){stop("Scaling must be 1 or a multiple of 3 (check MASSPRF file for details)")}
  MASSPRF_FILE<-read.csv(MASSPRF_Table,sep="\t")
  structUnaligned<-paste0(pdbseq(read.pdb(pdbFile)),collapse = "")
  origUnaligned<-gsub("\\*","",as.character(Biostrings::translate(Biostrings::readDNAStringSet(MASSPRF_Nuc_Fasta))[[1]]))
  tmpAln<-c(">orig",origUnaligned,">struct",structUnaligned)
  writeLines(tmpAln,con="_tmpAln.fasta")
  res<-muscle(Biostrings::readAAStringSet("_tmpAln.fasta"))
  file.remove("_tmpAln.fasta")
  res<-as.character(res)
  origSeq<-res[1]
  structSeq<-res[2]
  scaling<-scaling
  onlySig<-onlySig
  sigSetting<-sigSetting
  
  redHex<-rgb(rgb1[1],rgb1[2],rgb1[3],maxColorValue = 255)
  blueHex<-rgb(rgb2[1],rgb2[2],rgb2[3],maxColorValue = 255)
  numberOfColors<-bins
  ehHex<-rgb(ehColor[1],ehColor[2],ehColor[1],maxColorValue = 255)
  
  isSignificant<-function(LCI,UCI){
    if(is.na(LCI)){return(FALSE)}
    if(is.na(UCI)){return(FALSE)}
    if(LCI>0){return(TRUE)}
    if(UCI<0){return(TRUE)}
    return(FALSE)
  }
  
  
  
  
  tmpMASS<-MASSPRF_FILE[,c("Position","Gamma","Lower_CI_Gamma","Upper_CI_Gamma")]
  colnames(tmpMASS)<-c("Position","Gamma","LCI","UCI")
  GammaValues <- tmpMASS[, c("Position", "Gamma")]
  values<-GammaValues$Gamma
  #stolen from stackOverflow
  ii <- cut(values, breaks = seq(min(values), max(values), len = numberOfColors), 
            include.lowest = TRUE)
  ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
  myColors <- colorRampPalette(c(blueHex,redHex))((numberOfColors-1))[ii]
  tmpMASS[,"color"]<-myColors
  Raw<-data.frame(matrix(nrow=nchar(origUnaligned),ncol=6))
  colnames(Raw)<-c("AA","Gamma","LCI","UCI","color","Significant")
  Raw$AA<-unlist(strsplit(x = origUnaligned,split=""))
  
  i<-1
  if(scaling==1){
    meanAAPosGamma<-data.frame(matrix(nrow=nchar(origUnaligned),ncol=6))
    colnames(meanAAPosGamma)<-c("Position","Gamma","LCI","UCI","color","Significant")
    meanAAPosGamma$Position<-((1:nchar(origUnaligned))-1)
    massPRFCOPY<-tmpMASS[,c("Gamma","LCI","UCI","color")]
    for(j in 1:nrow(meanAAPosGamma)){
      theseThree<-massPRFCOPY[1:3,]
      theseThree<-theseThree[complete.cases(theseThree),]
      meanGam<-mean(theseThree$Gamma)
      meanLCI<-mean(theseThree$LCI)
      meanUCI<-mean(theseThree$UCI)
      meanSig<-NA
      sigs<-mapply(theseThree$LCI,FUN=isSignificant,UCI=theseThree$UCI)
      if(sigSetting=="any"){
        if(sum(sigs)>=1){meanSig<-TRUE}else{meanSig<-FALSE}
      }
      if(sigSetting=="majority"){
        if(sum(sigs)>=2){meanSig<-TRUE}else{meanSig<-FALSE}
      }
      if(sigSetting=="strict"){
        if(FALSE%in%sigs){meanSig<-FALSE}else{meanSig<-TRUE}
      }
      if(sigSetting=="average"){
        meanSig<-isSignificant(meanLCI,meanUCI)
      }
      meanAAPosGamma[j,]<-c(meanAAPosGamma$Position[j],meanGam,meanLCI,meanUCI,NA,meanSig)
      massPRFCOPY<-massPRFCOPY[-c(1,2,3),]
    }
    Raw$Gamma<-meanAAPosGamma$Gamma
    Raw$LCI<-meanAAPosGamma$LCI
    Raw$UCI<-meanAAPosGamma$UCI
    Raw$Significant<-meanAAPosGamma$Significant
    
    values<-Raw$Gamma
    #stolen from stackOverflow
    ii <- cut(values, breaks = seq(min(values), max(values), len = numberOfColors), 
              include.lowest = TRUE)
    myColors <- colorRampPalette(c(blueHex,redHex))((numberOfColors-1))[ii]
    Raw$color<-myColors
  }else{
    if(scaling%%3!=0){stop("needs to be divisible by 3!!!!")}
    aaScale<-scaling/3
    while(i < nrow(tmpMASS)) {
      thisG <- tmpMASS$Gamma[i]
      thisC <- tmpMASS$color[i]
      thisL<-tmpMASS$LCI[i]
      thisU<-tmpMASS$UCI[i]
      thisSig<-isSignificant(thisL,thisU)
      relevantIndecies<-((((i*aaScale)-(aaScale-1))):(i*aaScale))
      Raw$Gamma[relevantIndecies]<-thisG
      Raw$color[relevantIndecies]<-thisC
      Raw$LCI[relevantIndecies]<-thisL
      Raw$UCI[relevantIndecies]<-thisU
      Raw$Significant[relevantIndecies]<-thisSig
      i <- i + 1
    }
  }
  Align<-data.frame(matrix(ncol=5,nrow=nchar(structSeq)))
  colnames(Align)<-c("orig","struct","gamma","color","Significant")
  
  Align$orig<-unlist(strsplit(origSeq,""))
  Align$struct<-unlist(strsplit(structSeq,""))
  
  
  aIndex<-1
  rIndex<-1
  while(rIndex<nrow(Raw)){
    if(Align$orig[aIndex]==Raw$AA[rIndex]){
      Align$gamma[aIndex]<-Raw$Gamma[rIndex]
      Align$color[aIndex]<-Raw$color[rIndex]
      Align$Significant[aIndex]<-Raw$Significant[rIndex]
      rIndex<-rIndex+1
    } else if(Align$orig[aIndex]=="-"){
      Align$gamma[aIndex]<-NA
      Align$color[aIndex]<-ehHex #grey for missing
    }
    aIndex<-aIndex+1
  }
  
  fStruct<-Align[(Align$struct!="-"),]
  fStruct[is.na(fStruct$gamma),"color"]<-ehHex
  fStruct[is.na(fStruct$Significant),"Significant"]<-FALSE
  if(onlySig){
    fStruct[which(fStruct$Significant==FALSE),"color"]<-ehHex
  }
  
  
  setmyColor<-function(HEX,COLORNAME,POS){
    myString<- "colordef COLORNAME HEX; color COLORNAME #0:POS;"
    myString<-gsub("HEX",HEX,myString)
    myString<-gsub("COLORNAME",COLORNAME,myString)
    myString<-gsub("POS",POS,myString)
    return(myString)
  }
  
  holdAllCommands<-""
  for(i in 1:nrow(fStruct)){
    myHex<-fStruct$color[i]
    myPOS<-(i-1) #convert to python indexing
    tmp<-setmyColor(myHex,paste0("Custom",(myPOS)),myPOS)
    holdAllCommands<-c(holdAllCommands,tmp)
  }
  holdAllCommands<-holdAllCommands[nchar(holdAllCommands)>1]
  writeLines(holdAllCommands,con=outfile)
  
  
}



