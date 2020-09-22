#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#source("https://bioconductor.org/biocLite.R")  
#biocLite("Biostrings")

library(shiny)
library(shinydashboard)
library(shinyBS)
library(markdown)
library(seqinr)
library(Biostrings)
library(ShortRead)
library(stringr)
library(data.table)
library(DT)

# Define UI for application that draws a histogram
ui <- dashboardPage(skin="green",
                    
                    # Application title
                    # titlePanel("Frietze Lab Genomics"),
                    
dashboardHeader(title = 'Frietze Lab Genomics'),
   
  dashboardSidebar(
   # position="center",
  width=350,
  sidebarMenu(
      #column(1,

   fileInput("file1", "Choose Reference File (FASTA)" #,
   #accept = c(
  #   "text/csv",
    # "text/comma-separated-values,text/plain",
    # ".csv")
       ),

   #EcoRI: GAATTC
   #HindIII: AAGCTT
   #BglII: AGATCT
   #DpnII: GATC
   #Csp6I: GTAC
   #NlaIII: CATG
      div(
      selectInput("enzyme1", "Restriction Enzyme #1:",
               c("EcoRI - GAATTC" = "GAATTC",
                 "HindIII - AAGCTT" = "AAGCTT",
                 "BglII - AGATCT" = "AGATCT",
                 "DpnII - GATC" = "GATC",
                 "Csp6I - GTAC" = "GTAC",  
                 "NlaIII - CATG" = "CATG"
               )),
      style="margin-top:-30px;"),
   
   div(
   selectInput("enzyme2", "Restriction Enzyme #2:",
               c("EcoRI - GAATTC" = "GAATTC",
                 "HindIII - AAGCTT" = "AAGCTT",
                 "BglII - AGATCT" = "AGATCT",
                 "DpnII - GATC" = "GATC",
                 "Csp6I - GTAC" = "GTAC",  
                 "NlaIII - CATG" = "CATG"
               )),
   style="margin-top:-30px;"),
   
   div(
   sliderInput("fragment_1_size", "Fragment 1 size:",
               min = 500, max = 1500,
               value = 500),
   style="margin-top:-30px;"
   ),
   
   div(
   sliderInput("fragment_2_size", "Fragment 2 size:",
               min = 300, max =500 ,
               value = 300),
   style="margin-top:-30px;"
   ),
   
   #tags$head(tags$script(src = "message-handler.js")),
   actionButton("do", 
                "Get Primer Sequences",
                style="width:300px; 
                height:50px;
                font-size:20px;
                margin-top:0px;")
   ,
   downloadButton('downloadData', 
                'Download Primer Sequences', 
                style="width:300px; 
                height:50px; 
                font-size:20px; 
                text-align:center;
                margin-left:16px;
                margin-top:-50px;")
    ,
    downloadButton('downloadBed', 
               'Download Viewpoints Bed File', 
               style="width:300px; 
                height:50px; 
                font-size:20px; 
                text-align:center;
                margin-left:-300px;
               margin-top:100px;")
      )),

dashboardBody(

  fluidRow(align="center" ,
           
           ##tags$h1(
           ## style="font-size:16px;margin:40px;",
           ## "Welcome to ",tags$b("4C-PRIMER"), 
           ## "an app to create sequences that can
           ##  be used in conjunction with primer3 to design primers for 4C-sequencing
           ##  experiments.  Use the sidebar to fill out parameters and get the primer sequences.
           ##  The primer3 input sequences in the last column of the output (output can be downloaded 
           ##  by clicking the download button!) can then be entered 
           ##  into primer3 (", tags$a(href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/",
           ##                           "https://www.ncbi.nlm.nih.gov/tools/primer-blast/"), 
           ##  ") to obtain the reverse primers. Note that the first 20 bases of the primer3 input 
           ##  sequences are the forward primers."),
           
           tags$h1(
              style="font-size:20px;margin:40px;width:600px;",
              "Welcome to ",tags$b("4C-PRIMER"), 
              "an app to create sequences that can
              be used in conjunction with primer3 (", 
              tags$a(href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/",
              "https://www.ncbi.nlm.nih.gov/tools/primer-blast/"),
              ") to design primers for 4C-sequencing experiments."),
              
          # #htmlOutput(outputId = "image1",height=300, width=600)#,
          # #htmlOutput(outputId = "image0"),
          # #htmlOutput(outputId = image1")
          # #includeHTML("C:\\Users\\Mike\\Desktop\\shiny_apps\\cmb_retreat\\www\\hhv3_cohrs_multiqc_report.html")
          # #tableOutput("contents")
          # 
          # #htmlOutput("contents") 
          # 
          textOutput("output_seq_name"),
          textOutput("output_seq_length"),
          br(),
          DT::dataTableOutput(
            "contents",
            width="600px",
            height="300px"
          )
  )
  
))
      

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
########################### Declare helper functions #########################################
   
   dna_complement <-function(x){
     count<-1
     compl<-character(length=nchar(x))
     dna_string_split <- unlist(strsplit(x,""))
     for(i in dna_string_split){
       if(i=="G"){
         compl[count]<-"C"
       }else if(i=="C"){
         compl[count]<-"G"
       }else if(i=="A"){
         compl[count]<-"T"
       }else if(i=="T"){
         compl[count]<-"A"
       }else{
         stop("unknown nucleotide")
       }
       count=count+1
     }
     ##print(compl)
     return(paste0(compl,collapse=""))
   }
   
#dna_complement("ATCGCTAGCATGGCAAGTGTAAGAAAAAAAAAACCACGTTCCAAAAACACCTATACGGGAAACATCTCTTCACTCCCTCCCCGTCCCAACCACCGCAACACGAACTACGGTAGAAGCTTATCGCTAGCATGGCAAGTGTAAGAAAAAAAAAACCACGTTCCAAAAACACCTATACGGGAAACATCTCTTCACTCCCTCCCCGTCCCAACCACCGC")

########################### Run Script #################################################

rv <- reactiveValues(
  m=data.frame(),
  bed=data.frame()
)

observeEvent(input$do, {
  
   withProgress(message = 'Program has begun running...', value = 0, {
    
   inFile <- input$file1
   
   ##print(inFile$datapath)
  
   ##if(is.null(inFile))
   ##{
   ##  return(NULL)
   ##}
  
   ## print("testing ...")
   ## Read FASTA (or FASTQ) files in an XStringSet object:
   ##reads_in <- Biostrings::readFasta (filepath=inFile$datapath,format="fasta")
   reads_in <- Biostrings::readDNAStringSet(filepath=inFile$datapath,format="fasta")
   ##print(typeof(reads_in))
   
   ##print("tesint...")
   
   incProgress(1/5, detail = "Finished Loading Fasta File")
    
   ##library("Biostrings")   
   ##colnames(seq_data)<-names(reads_in)
   
   seq_string <- as.character(reads_in)
   seq_string <- tolower(stringr::str_squish(seq_string))
   
   seq_length<-nchar(unname(seq_string))
   
   output$output_seq_name<-renderText(
   {
      paste0("Seq name: ", names(reads_in)) 
   })
   
   output$output_seq_length<-renderText(
   {
      paste0("Seq length: ", seq_length)
   })
   
   incProgress(2/5, detail = "Gathered Chromosomes")
   
   ##Now let's get restriction sites:
   
   ##EcoRI: GAATTC
   ##eco_split<-unlist(strsplit(gsub("GAATTC"," ",seq_string),split=" "))
   
   ##HindIII: AAGCTT
   ##it cuts between the 2 A's. 
   ##print(input$enzyme1)
   ##print(input$enzyme2)
   
   enz_matches <- str_locate_all(seq_string, tolower(input$enzyme1))
   ##Get all start positions (which are the same as end positions) 
   ##of the "_" chars
   enz_splits <- enz_matches[[1]]
   
   frag_starts<-numeric(0)
   frag_stops<-numeric(0)
   frag_seq<-character(0)
   frag_length<-numeric(0)
   
   for(i in seq(from=1,by=1,to=nrow(enz_splits)))
   {
     if(i==1)
     {
       frag_start <- 1
       frag_stop <- enz_splits[1,2]
       frag_seq <- substr(seq_string, frag_start, frag_stop)
       frag_length <- nchar(frag_seq)
     }else {
       new_start <- enz_splits[i-1,2]+1
       new_stop <- enz_splits[i,2] 
       new_seq <- substr(seq_string, new_start, new_stop)
       new_length <- nchar(new_seq)
       frag_start <- c(frag_start, new_start)
       frag_stop <- c(frag_stop, new_stop) 
       frag_seq <- c(frag_seq, new_seq)
       frag_length <- c(frag_length, new_length)
     }
    ##print(length(enz_matches))
   }
   
   ##hind_split<-unlist(strsplit(seq_string,split=tolower(input$enzyme1)))
   
   first_frame<-data.frame(
      seq=character(length(frag_length)),
      start=numeric(length(frag_length)),
      stop=numeric(length(frag_length)),
      length=numeric(length(frag_length)),
      stringsAsFactors = FALSE
      )
   
   first_frame$seq <- frag_seq
   first_frame$start <- frag_start
   first_frame$stop <- frag_stop
   first_frame$length <- frag_length
   
   if(!grepl("aagctt",first_frame[nrow(first_frame),1]))
   {
     ##print("true!")
     first_frame<-first_frame[-nrow(first_frame),]
   }
   
   incProgress(3/5, detail = "Finished First Enzyme Cut")
   
   ##hind_split<-hind_split[nchar(hind_split)>=input$viewpoint_size]
   first_frame<-first_frame[first_frame$length>=input$fragment_1_size,]
   
   second_frame <- data.frame(
     first_seq = character(0),
     first_start = numeric(0),
     first_stop = numeric(0),
     first_length = numeric(0),
     second_seq = character(0),
     second_start = numeric(0),
     second_stop = numeric(0),
     second_length =numeric(0),
     forward_primer = character(0),
     sequence_for_primer3 = character(0),
     stringsAsFactors = FALSE
     )
   
   for(j in seq(from=1,by=1,to=nrow(first_frame)))
   {
     ##print(seq(from=1,by=1,to=nrow(first_frame)))
     enz2_matches <- str_locate_all(tolower(first_frame[j,1]), 
                                    tolower(input$enzyme2))
     enz2_splits <- enz2_matches[[1]]
     ##If there are no second enzyme sites, skip to next:
     if(length(enz2_splits)==0){
       next
     }
     
     ##Get the right-most fragment:
     enz2_split <- enz2_splits[nrow(enz2_splits),]
     
     frag_seq_2 <- substr(first_frame[j,1], 
                          enz2_split[1], 
                          nchar(first_frame[j,1]))
     frag_length_2 <- nchar(frag_seq_2)
     frag_stop_2 <- first_frame[j,3]
     frag_start_2 <- frag_stop_2 - frag_length_2 + 1  
     
     second_int_frame<-data.frame(
     first_seq = first_frame[j,1],
     first_start = first_frame[j,2],
     first_stop = first_frame[j,3],
     first_length = first_frame[j,4],
     second_seq = frag_seq_2,
     second_start = frag_start_2,
     second_stop = frag_stop_2,
     second_length = frag_length_2,
     forward_primer = "",
     sequence_for_primer3 = "",
     stringsAsFactors = FALSE)
     
     second_frame<-rbind(second_frame,second_int_frame)
   }
     
  ##Filter second digest fragments by minimum size:
     
  second_frame<-second_frame[second_frame$second_length>=input$fragment_2_size,]
     
  incProgress(4/5, detail = "Finished Second Enzyme Cut")

  ##Output the table
  second_frame$sequence_for_primer3<-paste0(
     substr(second_frame$second_seq,
            nchar(second_frame$second_seq)-19,
            nchar(second_frame$second_seq)),
     paste(rep("N",times=100),collapse=""),
     substr(second_frame$second_seq,1,100)
   )
   
   ##frag_frame$forward_primer<-sapply(X=substr(frag_frame$primer_seq,1,20),FUN=dna_complement)
   second_frame$forward_primer<-substr(second_frame$sequence_for_primer3,1,20)
   
   ##Can output to console for debugging:
   ##tab-delimited:
   ##cat(paste(second_frame$sequence_for_primer3,
   ##          second_frame$forward_primer,
   ##          second_frame$second_seq,
   ##          second_frame$second_length,
   ##          second_frame$second_start,
   ##          second_frame$second_stop,
   ##          sep="\t"), 
   ##          sep="\n")
   
   ##cat(frag_frame$viewpoint_fragment_length,sep="\n")
   ##For reverse primer enter each primer_seq into primer3. 
   rv$m<-second_frame
   
   bed<-second_frame[,c(6,7)]
   bed$chrom<-rep(names(reads_in),times=nrow(second_frame))
   colnames(bed)<-c("chromStart","chromEnd","chrom")
   bed$chromStart<-bed$chromStart - 1
   bed$chromEnd<-bed$chromEnd
   bed$length <- bed$chromEnd - bed$chromStart 
   bed<-bed[,c("chrom","chromStart","chromEnd","length")]
   rv$bed<-bed
   
   screen_frame <- second_frame[,c("forward_primer","sequence_for_primer3")]
   rv$screen_frame <- screen_frame
   
   output$contents <- DT::renderDataTable(
     {rv$screen_frame},
     options = list(scrollX = TRUE, scrollY=TRUE),
     rownames= FALSE
     ##thedata(frag_frame)
     ##frag_frame$primer_seq
     )

   incProgress(5/5, detail = "Finished Output!") 
   
   })
   
})

output$downloadData <- downloadHandler(
  filename = function() {
    ##paste(input$dataset, ".csv", sep = "")
    paste("sequences_for_primer3", ".csv", sep = "")
  },
  content = function(file) {
    ##Can output more comprehensive output:
     ##write.csv(rv$m, file, row.names = FALSE)
     
     ##Or just output the forward and reverse primers:
     write.csv(rv$screen_frame, file, row.names = FALSE)
   }
 )

output$downloadBed <- downloadHandler(
  filename = function() {
    ##paste(input$dataset, ".csv", sep = "")
    paste("viewpoints", ".bed", sep = "")
  },
  content = function(file) {
    write.table(rv$bed, file, row.names = FALSE,col.names=FALSE, quote=FALSE)
  }
)

}

# Run the application 
shinyApp(ui = ui, server = server)
