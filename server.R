shinyServer(function(input, output) {
  
  DATA1 <- reactive({
    if (is.null(input$file1)) {return()} else {
      dat <- read.fasta(input$file1$datapath) #路徑為input$file1$datapath
      return(dat) 
    }
  })
  
  DATA2 <- reactive({
    if (is.null(input$file2)) {return()} else {
      dat <- read.fasta(input$file2$datapath) #路徑為input$file2$datapath
      return(dat) 
    }
  })
  
  ALIGNOUTCOME <- eventReactive(input$submit,{
    seq1_list = DATA1()
    seq2_list = DATA2()
    if (is.null(seq1_list)|is.null(seq2_list)) {return()} else {
      seq1_seqs=getSequence(seq1_list) 
      seq2_seqs=getSequence(seq2_list) 
      
      seq1_names=getName(seq1_list) 
      seq2_names=getName(seq2_list) 
      
      seq1_annot=getAnnot(seq1_list)
      seq2_annot=getAnnot(seq2_list)
      
      seq1_seqs_unique=seq1_seqs[!duplicated(seq1_seqs)] 
      seq2_seqs_unique=seq2_seqs[!duplicated(seq2_seqs)]
      
      seq1_names_unique=seq1_names[!duplicated(seq1_seqs)]
      seq2_names_unique=seq2_names[!duplicated(seq2_seqs)]
      
      seq1_annot_unique=seq1_annot[!duplicated(seq1_seqs)]
      seq2_annot_unique=seq2_annot[!duplicated(seq2_seqs)]
      
      seq1_list_unique=list() 
      
      for(a in 1:length(seq1_seqs_unique)){
        mseq=seq1_seqs_unique[[a]]
        mname=seq1_names_unique[[a]]
        mannot=seq1_annot_unique[[a]]
        seq1_list_unique[[a]]=as.SeqFastadna(mseq, name =mname , Annot =mannot) 
      }
      
      seq2_list_unique=list()
      for(b in 1:length(seq2_seqs_unique)){
        mseq=seq2_seqs_unique[[b]]
        mname=seq2_names_unique[[b]]
        mannot=seq2_annot_unique[[b]]
        seq2_list_unique[[b]]=as.SeqFastadna(mseq, name =mname , Annot =mannot) 
      }
      
      alignoutcome <- data.frame(seq1Name=character(),
                                 seq2Name=character(),
                                 seq1Annot=character(),
                                 seq2Annot=character(),
                                 seq1Len=double(),
                                 seq2Len=double(),
                                 alignscore=double(),
                                 possibleMaxScore=double(),
                                 stringsAsFactors=FALSE)
      
      data(BLOSUM50) 
      
      idx=1
      aligntimes=printAlignNo(seq1_list,seq2_list)
      n=aligntimes[[3]]
      
      withProgress(message = "In processing...",value=0,{
        for (i in 1:length(seq1_list_unique)){
          seq1_element=seq1_list_unique[[i]]
          seq1_string=convert_string_upper(seq1_element) #將蛋白質序列轉換成英文大寫及字串形式，方能進行後續的序列比對
          
          for (j in 1:length(seq2_list_unique)){
            seq2_element=seq2_list_unique[[j]]
            seq2_string=convert_string_upper(seq2_element) #將蛋白質序列轉換成英文大寫及字串形式，方能進行後續的序列比對
            
            alignmentscore =pairwiseAlignment(seq1_string,seq2_string,substitutionMatrix=BLOSUM50,gapOpening=-2,gapExtension=-8,scoreOnly=TRUE,type="local")
            alignoutcome[idx,"seq1Name"]= getName(seq1_element) #目前被比對之seq1序列名稱
            alignoutcome[idx,"seq2Name"]= getName(seq2_element) #目前被比對之seq2序列名稱
            alignoutcome[idx,"seq1Annot"]= getAnnot(seq1_element) #目前被比對之seq1序列註解
            alignoutcome[idx,"seq2Annot"]= getAnnot(seq2_element) #目前被比對之seq2序列註解
            alignoutcome[idx,"seq1Len"]= getLength(seq1_element) #seq1序列長度
            alignoutcome[idx,"seq2Len"]= getLength(seq2_element) #seq2序列長度
            alignoutcome[idx,"alignscore"]= alignmentscore #目前被比對之score   
            alignoutcome[idx,"possibleMaxScore"]= alignPossibleMaxScore(seq1_string,seq2_string)#目前被比對之score          
            idx=idx+1
            incProgress(1/n)
          }
        }
      })
      
      alignoutcome$matchscale=alignoutcome$alignscore/alignoutcome$possibleMaxScore
      return(alignoutcome)
    }
  })
  
  output$choose_columns1 <- renderUI({
    alignoutcome = ALIGNOUTCOME()
    if (is.null(alignoutcome)) {return()} else {
      MIN = min(alignoutcome$alignscore)
      MAX = max(alignoutcome$alignscore)
      VALUE = quantile(alignoutcome$alignscore, 0.25)
      sliderInput("cut", label = "", min = MIN, max = MAX, value = VALUE)
    }
  })
  
  PIECHARTDATA <- reactive({
    alignoutcome = ALIGNOUTCOME()
    if (is.null(alignoutcome)|is.null(input$cut)) {return()} else {
      dat = PieChartData(alignoutcome,input$cut)
      return(dat)
    }
  })
  
  output$plot1 <- renderGvis({
    PieData = PIECHARTDATA()
    if (is.null(PieData)) {return()} else {
      Pie = gvisPieChart(PieData[[1]], options=list(width=input$MonitorWidth/2-150, height=input$MonitorHeight-400))
      Pie
    }
  })
  
  output$plot2 <- renderGvis({
    PieData = PIECHARTDATA()
    if (is.null(PieData)) {return()} else {
      Pie = gvisPieChart(PieData[[2]], options=list(width=input$MonitorWidth/2-150, height=input$MonitorHeight-400))
      Pie
    }
  })
  
  output$download1 <- downloadHandler(
    filename = function() {'Result.txt'},
    content = function(con) {
      alignoutcome = ALIGNOUTCOME()
      if (is.null(alignoutcome)|is.null(input$cut)) {return()} else {
        write.table(alignoutcome[alignoutcome$alignscore>input$cut,], con, sep="\t", row.names=FALSE, col.names=TRUE)
      }
    }
  )
  
  output$download2 <- downloadHandler(
    filename = function() {'plot.pdf'},
    content = function(con) {
      PieData = PIECHARTDATA()
      if (is.null(PieData)) {return()} else {
        pdf(con, family = "serif")
        for (j in 1:length(PieData)) {
          pie.scales = PieData[[j]][,2]
          names(pie.scales) = PieData[[j]][,1]
          pie(pie.scales, col = rainbow(length(pie.scales)))
        }
        dev.off()
      }
    }
  )

  
})