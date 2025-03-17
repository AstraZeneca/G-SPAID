# 
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFiles)
library(spatstat)
library(DescTools)
library(dplyr)
library(reshape2)
library(rstatix)
library(purrr)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(circlize)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("G-SPAID"),
  
  fileInput("file1", "Choose CSV file", multiple = FALSE,
            accept = c("text/csv","text/comma-separated-values,text/plain",".csv"),
            width = NULL, buttonLabel = "Browse...",
            placeholder = "No file selected"),
  
  selectInput("series", "Cell type:", choices=c(),multiple=TRUE),
  
  sidebarPanel(
    selectInput("Imagecol", "Image", choices=c(),multiple=FALSE),
    selectInput("group", "group Column", choices=c(),multiple=FALSE),
    
    selectInput("Xmin", "X min Column", choices=c(),multiple=FALSE),
    selectInput("Xmax", "X max Column", choices=c(),multiple=FALSE),
    selectInput("Ymin", "Y min Column", choices=c(),multiple=FALSE),
    selectInput("Ymax", "Y max Column", choices=c(),multiple=FALSE),
    
    numericInput('sim','no sim',20),
    numericInput('r','distance',30),
    numericInput('m','rank',5),
    numericInput('ncell','min cell count',10),
    
    
    actionButton("go","run")
  ),
  
  mainPanel(
    textOutput('dir'),
    #InteractiveComplexHeatmapOutput("ht1")
    
    
    
  )
)

#################################################################################
#
# Server function
#
server <- function(input, output) {
  options(shiny.maxRequestSize=250*1024^2)
  
  # path
  path <- reactiveVal(getwd())
  
  input_df <- reactive({
    req(input$file1)
    data.frame(read.csv(input$file1$datapath))    
  })
  
  observeEvent(input$file1,{
    updateSelectInput(inputId = "series",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "Imagecol",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "group",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "Xmin",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "Xmax",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "Ymin",
                      choices = colnames(input_df()))
  })
  observeEvent(input$file1,{
    updateSelectInput(inputId = "Ymax",
                      choices = colnames(input_df()))
  })
  
  
  observeEvent(input$go,{
    
    df <- input_df()[,c(input$series,input$Imagecol,input$Xmin,input$Xmax,input$Ymin,input$Ymax,input$group)]
    colnames(df) = c(input$series,'Image','Xmin','Xmax','Ymin','Ymax','Group')
    
    if (is.numeric(df$Xmin) & is.numeric(df$Xmax) & is.numeric(df$Ymin) & is.numeric(df$Ymax) == TRUE){
      df$X <- (df$Xmin + df$Xmax)/2
      df$Y <- (df$Ymin + df$Ymax)/2
      dt <- df[,c(input$series,"Image","X","Y","Group")]
      
      path(file.path(path(),
                     paste0('interaction_',
                            gsub(':', '', format(Sys.time(), "%d%m%y%X")))))
      dir.create(path())
      outDir <- path()
      
      #normfname <- file.path(outDir,("input_data.csv"))
      #write.csv(dt, normfname,row.names=FALSE)
      
      print("Calculating ...")
      print({paste("Image no.", length(unique(input_df()$Image)), ",", "no. cell type", length(input$series), ",", "distance range", input$distance)})
      
      calculate_AUC <- function(index, d) {
        
        gd <- g$obs[0:index]
        gd_hi <- g$hi[0:index]
        gd_lo <- g$lo[0:index]
        g30 <- g %>% filter(r <= input$r)
        
        if (isTRUE(all(gd == gd_hi|gd == gd_lo ))) {
          #print("all gd not sig diff")
          AUC_cluster <<- c(AUC_cluster,0)
          AUC_dispersion <<- c(AUC_dispersion,0)
          dist_c <<- c(dist_c,paste(0,d,"insig"))
          dist_d <<- c(dist_d,paste(0,d,"insig"))
          
        } else if (isTRUE(all(gd >= gd_hi))) {
          #print("all gd greater that gd_hi")
          auc <- AUC(g$r,g$obs, to = d)/(d)
          theo <- AUC(g$r,g$hi, to = d)/(d)
          AUC_cluster <<- c(AUC_cluster,auc-theo)
          AUC_dispersion <<- c(AUC_dispersion,0)
          
          dist_c <<- c(dist_c,paste(0,d,"gd>hi"))
          dist_d <<- c(dist_d,paste(0,d,"gd>hi"))
          
        } else if (isTRUE(all(gd <= gd_lo))) {
          #print("all gd less that gd_lo")
          auc <- AUC(g$r,g$obs, to = d)/(d)
          theo <- AUC(g$r,g$lo, to = d)/(d)
          
          AUC_cluster <<- c(AUC_cluster,0)
          AUC_dispersion <<- c(AUC_dispersion,theo-auc)
          
          dist_c <<- c(dist_c,paste(0,d,"gd<lo"))
          dist_d <<- c(dist_d,paste(0,d,"gd<lo"))
          
        } else if (isTRUE(all(gd <= gd_hi & gd >= gd_lo))) {
          #print("all gd not sig diff")
          AUC_cluster <<- c(AUC_cluster,0)
          AUC_dispersion <<- c(AUC_dispersion,0)
          dist_c <<- c(dist_c,paste(0,d,"insig"))
          dist_d <<- c(dist_d,paste(0,d,"insig"))
          
        } else if (isTRUE(any(gd > gd_hi)) & (any(gd < gd_lo))) {
          #print("gd greater&less than theo at some r range")
          g30$logiclo <- g30$lo - g30$obs
          runs <- rle(g30$logiclo > 0)          
          Trun <- which(runs$values == TRUE)    
          runs.lengths.cumsum = cumsum(runs$lengths)
          ends = runs.lengths.cumsum[Trun]
          newindex = ifelse(Trun>1, Trun-1, 0)
          starts = runs.lengths.cumsum[newindex] + 1
          if (0 %in% newindex) starts = c(1,starts)
          
          difflo <- c()
          for (runi in 1:length(starts)){
            dlo_from <- g30$r[starts[runi]]
            dlo_to <- g30$r[ends[runi]]
            
            auc <- AUC(g30$r,g30$obs,from = dlo_from, to = dlo_to)/(d)
            theo <- AUC(g30$r,g30$lo,from = dlo_from, to = dlo_to)/(d)
            difflo <- c(difflo,theo - auc)
          }
          
          g30$logicho <- g30$obs - g30$hi
          runs <- rle(g30$logicho > 0)         
          Trun <- which(runs$values == TRUE) 
          runs.lengths.cumsum = cumsum(runs$lengths)
          ends = runs.lengths.cumsum[Trun]
          newindex = ifelse(Trun>1, Trun-1, 0)
          starts = runs.lengths.cumsum[newindex] + 1
          if (0 %in% newindex) starts = c(1,starts)
          
          diffhi <- c()
          for (runi in 1:length(starts)){
            dlo_from <- g30$r[starts[runi]]
            dlo_to <- g30$r[ends[runi]]
            
            auc <- AUC(g30$r,g30$obs,from = dlo_from, to = dlo_to)/(d)
            theo <- AUC(g30$r,g30$lo,from = dlo_from, to = dlo_to)/(d)
            diffhi <- c(diffhi,abs(auc -theo))
          }
          
          AUC_cluster <<- c(AUC_cluster,sum(diffhi))
          AUC_dispersion <<- c(AUC_dispersion,sum(difflo))
          dist_c <<- c(dist_c,"both")
          dist_d <<- c(dist_d,"both")
          
        } else if (isTRUE(any(gd < gd_lo))) {
          #print("gd less than theo at some r range")
          g30$logic <- g30$lo - g30$obs
          runs <- rle(g30$logic > 0)         
          Trun <- which(runs$values == TRUE) 
          runs.lengths.cumsum = cumsum(runs$lengths)
          ends = runs.lengths.cumsum[Trun]
          newindex = ifelse(Trun>1, Trun-1, 0)
          starts = runs.lengths.cumsum[newindex] + 1
          if (0 %in% newindex) starts = c(1,starts)
          
          difflo <- c()
          dlo_froml <- c()
          dlo_tol <- c()
          runi <- 1
          for (runi in 1:length(starts)){
            dlo_from <- g30$r[starts[runi]]
            dlo_to <- g30$r[ends[runi]]
            
            dlo_froml <- c(dlo_froml,dlo_from)
            dlo_tol <- c(dlo_tol,dlo_to)
            
            auc <- AUC(g30$r,g30$obs,from = dlo_from, to = dlo_to)/(d)
            theo <- AUC(g30$r,g30$lo,from = dlo_from, to = dlo_to)/(d)
            difflo <- c(difflo,abs(theo-auc))
          }
          
          AUC_cluster <<- c(AUC_cluster,0)
          AUC_dispersion <<- c(AUC_dispersion,sum(difflo))
          dist_c <<- c(dist_c,paste(0,d,"insig"))
          dist_d <<- c(dist_d,paste(paste(min(round(dlo_froml)),max(round(dlo_tol)),"multi dispersion")))
          
        } else if (isTRUE(any(gd > gd_hi))) {
          #print("gd greater than theo at some r range")
          g30$logic <- g30$obs - g30$hi
          runs <- rle(g30$logic > 0)         ### get Consecutive runs
          Trun <- which(runs$values == TRUE) ### Consecutive TRUE
          runs.lengths.cumsum = cumsum(runs$lengths)
          ends = runs.lengths.cumsum[Trun]
          newindex = ifelse(Trun>1, Trun-1, 0)
          starts = runs.lengths.cumsum[newindex] + 1
          if (0 %in% newindex) starts = c(1,starts)
          
          diffhi <- c()
          dlo_froml <- c()
          dlo_tol <- c()
          
          runi<-1
          for (runi in 1:length(starts)){
            dlo_from <- g30$r[starts[runi]]
            dlo_to <- g30$r[ends[runi]]
            
            dlo_froml <- c(dlo_froml,dlo_from)
            dlo_tol <- c(dlo_tol,dlo_to)
            
            auc <- AUC(g30$r,g30$obs,from = dlo_from, to = dlo_to)/(d)
            theo <- AUC(g30$r,g30$lo,from = dlo_from, to = dlo_to)/(d)
            diffhi <- c(diffhi,abs(auc -theo))
          }
          
          AUC_cluster <<- c(AUC_cluster,sum(diffhi))
          AUC_dispersion <<- c(AUC_dispersion,0)
          dist_d <<- c(dist_d,paste(0,d,"insig"))
          dist_c <<- c(dist_c,paste(paste(min(round(dlo_froml)),max(round(dlo_tol)),"multi cluster")))
          
        } else {
          #print("no value calculated")
          AUC_cluster <<- c(AUC_cluster,"error")
          AUC_dispersion <<- c(AUC_dispersion,"error")
          dist_d <<- c(dist_d,"error")
          dist_c <<- c(dist_c,"error")
          
        }
      }
      
      count_cell <- input$ncell
      cell_type <- c(input$series)
      distance <- input$r    # distance
      sim <- input$sim # 199        # number of simulation
      m <- input$m            # confidence interval ie 5th highest/lowest value from simulation
      sample_ID <- c()
      group_ID <- c()
      AUC_cluster <- c()
      AUC_dispersion <- c()
      dist_c <- c()
      dist_d <- c()
      count = 0
      
      for (ID in unique(dt$Image)){
        count = count + 1
        sample_ID <- c(sample_ID,ID)
        print(paste(count,"/",length(unique(dt$Image))))
        #print(ID)
        win <- filter(dt,Image == ID)
        group_ID <- c(group_ID,unique(win$Group))
        ID_window <- as.owin(ripras(win$X,win$Y))
        ixj <- c()
        
        for (i in cell_type){
          p1x <- win$X[win[i] > 0]
          p1y <- win$Y[win[i] > 0]
          p1 <- ppp(p1x, p1y, window=ID_window)
          
          for (j in cell_type){
            print(paste(i,j))
            ixj <- c(ixj,paste(i,'x',j))
            
            p2x <- win$X[win[j] > 0]
            p2y <- win$Y[win[j] > 0]
            p2 <- ppp(p2x, p2y, window=ID_window)
            
            if((length(p1x) <= count_cell |length(p2x) <= count_cell) == TRUE) {
              dist_d <- c(dist_d,distance)
              AUC_cluster <- c(AUC_cluster,"low-cell")
              dist_c <- c(dist_c,distance)
              AUC_dispersion <- c(AUC_dispersion,"low-cell")
              
            } else if((i==j) == TRUE){
              #print("Gest")
              g <- envelope(p1, Gest, nsim = sim, fix.n = TRUE, fix.marks = TRUE, savefuns = TRUE, nrank = m, verbose = FALSE)
              
              if (max(g$r) < distance){
                dvalue <- max(g$r)
              } else {
                dvalue <- distance
              }
              
              index_value <- which.min(abs(g$r - dvalue))
              calculate_AUC(index = index_value, d = dvalue)
              
            } else {
              #print("Gcross")
              pm = ppp(c(p1x, p2x),c(p1y,p2y),window= ID_window,marks=factor(c(rep(i,p1$n),rep(j,p2$n))))
              g <- envelope(pm, Gcross, i = i, j = j, nsim = sim, fix.n = TRUE, fix.marks = TRUE, savefuns = TRUE, nrank = m, verbose = FALSE) # Gcross
              
              if (max(g$r) < distance){
                dvalue <- max(g$r)
              } else {
                dvalue <- distance
              }
              
              index_value <- which.min(abs(g$r - dvalue))
              calculate_AUC(index = index_value, d = dvalue)
              
            }
          }
        }
      }
      AUC_df <- as.data.frame(matrix(unlist(AUC_dispersion), ncol = (length(cell_type)^2), byrow = TRUE,dimnames = list(sample_ID, ixj)))
      normfname <- file.path(outDir,("gdispersion_table.csv"))
      write.csv(AUC_df, normfname)
      
      AUC_df <- as.data.frame(matrix(unlist(AUC_cluster), ncol = (length(cell_type)^2), byrow = TRUE,dimnames = list(sample_ID, ixj)))
      AUC_df$sample_ID <- rownames(AUC_df)
      save_file1 <- melt(AUC_df,id.vars = "sample_ID")
      
      dist_df <- as.data.frame(matrix(unlist(dist_d), ncol = (length(cell_type)^2), byrow = TRUE,dimnames = list(sample_ID, ixj)))
      normfname <- file.path(outDir,("gdispersion_dist_table.csv"))
      write.csv(dist_df, normfname,row.names = FALSE)
      
      dist_df <- as.data.frame(matrix(unlist(dist_c), ncol = (length(cell_type)^2), byrow = TRUE,dimnames = list(sample_ID, ixj)))
      dist_df$sample_ID <- rownames(dist_df)
      save_file2 <- melt(dist_df,id.vars = "sample_ID")
      save_file <- left_join(save_file1,save_file2, by=c("sample_ID","variable"))
      colnames(save_file) <- c("sample_ID", "interaction","g area","distance")
      
      if (length(unique(dt$Group)) == 2){
        
        #cell_count
        cell_count <- c()
        for (ID in unique(dt$Image)){
          win <- filter(dt,Image == ID)
          for (i in cell_type){
            cell_count <- c(cell_count,length(win$X[win[i] > 0]))
          }  
        }    
        df_n <- as.data.frame(t(matrix(unlist(cell_count), ncol = length(unique(dt$Image)), nrow = length(cell_type))))
        colnames(df_n) = cell_type
        df_n$group = group_ID
        df_count <- melt(df_n,id.vars = "group")
        df_n$sample_ID = sample_ID
        #normfname <- file.path(outDir,("cellcount_table.csv"))
        #write.csv(df_n, normfname,row.names = FALSE)
        #print(paste0("cellcount_table saved ",outDir))
        df_count <- suppressMessages(df_count %>%
                                       group_by(group,variable) %>% 
                                       summarise(count = mean(value,na.rm = TRUE)) %>%
                                       as.data.frame() %>% 
                                       arrange(match(variable, rep(cell_type,each = 2))) %>% 
                                       pull(count))
        count <- round(as.numeric(df_count),0)
        
        #AUC matrix
        AUC_df <- as.data.frame(matrix(unlist(AUC_cluster), ncol = (length(cell_type)^2), byrow = TRUE,dimnames = list(sample_ID, ixj)))
        AUC_df <- replace(AUC_df, AUC_df=="low-cell", NA)
        AUC_df <- replace(AUC_df, AUC_df=="Inf", NA)
        AUC_df <- as.data.frame(sapply(AUC_df,as.numeric))
        AUC_df$group <- group_ID
        AUC_df <- melt(AUC_df,id.vars = "group")
        mean_df <- suppressMessages(AUC_df %>%
                                      group_by(group,variable) %>% 
                                      summarise(mean_AUC = mean(value,na.rm = TRUE),sd = sd(value,na.rm = TRUE)) %>%
                                      as.data.frame())
        
        write.csv(AUC_df, normfname,row.names = FALSE)
        
        print("XX")
        stat.test <- suppressMessages(
          AUC_df %>% 
            na.omit() %>% 
            group_by(variable) %>%
            nest() %>%  # Nest the data by variable
            mutate(
              test_result = map(data, ~{
                # Check if both groups have sufficient values
                if (length(unique(.x$group)) < 2 || min(table(.x$group)) < 2) {
                  # If not enough values, return a default tibble with zeros
                  tibble(
                    p = NA_real_,
                    p.adj = NA_real_,
                    p.adj.signif = "ns"
                  )
                } else {
                  # Perform the test if sufficient values
                  .x %>%
                    wilcox_test(value ~ group) %>%
                    adjust_pvalue(method = "BH") %>%
                    add_significance("p.adj")
                }
              })
            ) %>%
            unnest(test_result) %>%
            select(variable, p, p.adj, p.adj.signif)
        )
        print("XX")
        
        stat.test1 <- stat.test
        stat.test$group <- unique(dt$Group)[1]
        stat.test1$group <- unique(dt$Group)[2]
        stat.test1$p.adj.signif <- "ns"
        stat.final <- rbind(stat.test,stat.test1)
        stat.mean <- left_join(mean_df,stat.final,by=c("variable","group"))
        stat.mean <- stat.mean %>% arrange(match(variable, rep(ixj,each = 2)),match(group, rep(unique(dt$Group),length(cell_type))))
        
        # heatmap
        mean_matrix <- stat.mean %>% select(mean_AUC)
        mean_matrix <- t(matrix(unlist(mean_matrix), nrow = (length(cell_type)*2), ncol = length(cell_type)))
        mean_matrix <- scale(mean_matrix)
        colnames(mean_matrix) <- rep(unique(dt$Group),length(cell_type))
        rownames(mean_matrix) <- cell_type
        order_col <- cell_type[rep(1:length(cell_type), each = 2)]
        
        matrix_lab <- t(matrix(unlist(stat.mean$p.adj.signif),nrow = (length(cell_type)*2),ncol = length(cell_type)))
        matrix_lab[matrix_lab== "1"] <- " "
        
        # saving file
        normfname <- file.path(outDir,("stat_table.csv"))
        stat.mean <- stat.mean %>% select(-c(p.adj.signif))
        write.csv(stat.mean, normfname,row.names = FALSE)
        print(paste0("stat saved ",outDir))
        
        top_an <- HeatmapAnnotation(count = anno_barplot(count,add_numbers = TRUE))
        factorp <- (1.5/24)*(length(cell_type)*2)
        f1 <- colorRamp2(seq(min(mean_matrix, na.rm=T), max(mean_matrix, na.rm=T), length = 3), c("blue", "#EEEEEE", "red"))
        
        pdffile <- file.path(outDir,('heatmap.pdf'))
        pdf(pdffile, width = 11, height = 8)
        print(Heatmap(mean_matrix, column_split = factor(order_col,levels = unique(order_col)),na_col = "grey", col = f1,
                      column_names_gp = gpar(fontsize = 8),column_title_rot = 45,column_title_gp = gpar(fontsize = 10,fontface="bold"),row_names_gp = gpar(fontsize = 10,fontface="bold"),
                      show_column_dend = FALSE, show_row_dend = FALSE,cluster_columns = FALSE,cluster_rows = FALSE,cluster_column_slices = FALSE,name = "AUC",row_title = "i", column_names_rot = 90,
                      top_annotation = top_an, cell_fun = function(j, i, x, y, width, height, fill) {
                        if(matrix_lab[i, j] < 0.05)
                          grid.text(matrix_lab[i, j], x/factorp, y, gp = gpar(fontsize = 13))
                      }
        ))
        dev.off()
        print("heatmap saved")
        
        sig_interaction <- stat.test %>% filter(p.adj < 0.05) %>% pull(variable)
        count = 0 
        
        for (int in sig_interaction){
          count = count + 1
          print(paste(count,"/",length(sig_interaction),"ploting sig int"))  
          ### create a folder of the interaction 
          subDir <- gsub('[\t\n/+.-]', '', int)
          dir.create(file.path(outDir, subDir), showWarnings = FALSE)
          setwd(file.path(outDir, subDir))
          
          i <- strsplit(int," x ")[[1]][1]
          j <- strsplit(int," x ")[[1]][2]
          r <- input$r 
          
          sample_ID <- c()
          p1n <- c()
          p2n <- c()
          average_nnd <- c()
          group_ID <- c()
          AUC_cluster <- c()
          AUC_dispersion <- c()
          dist_c <- c()
          dist_d <- c()
          
          for (diet in unique(dt$Group)) {
            sub <- dt %>% filter(Group == diet)
            
            nnd_density <- c()
            subject <- c()
            
            pdf(paste0(diet,"_g-%03d.pdf"),width=11.4, height=6.5,onefile = FALSE)
            par(mfrow = c(4, 7),mgp=c(1,0.3,0),mar = c(2,1.5,1.5,1.5),oma = c(1, 1, 0, 0))
            
            for (x in unique(sub$Image)){
              group_ID <- c(group_ID,diet)
              win <- filter(dt,Image == x)
              
              #CREATING DATAFRAME 
              sample_ID <- c(sample_ID,x)
              p1x <- win$X[win[i] > 0]
              p1y <- win$Y[win[i] > 0]
              p2x <- win$X[win[j] > 0]
              p2y <- win$Y[win[j] > 0]
              p1n <- c(p1n,length(p1x))
              p2n <- c(p2n,length(p2x))
              
              window <- as.owin(ripras(win$X,win$Y))
              p1 <- ppp(p1x, p1y, window=window)
              p2 <- ppp(p2x, p2y, window=window)
              
              if((length(p1x) <= 1 |length(p2x) <= 1) == TRUE) {
                dist_d <- c(dist_d,distance)
                AUC_cluster <- c(AUC_cluster,"low-cell")
                dist_c <- c(dist_c,distance)
                AUC_dispersion <- c(AUC_dispersion,"low-cell")
                
                average_nnd <- c(average_nnd,"low-cell")
                nnd_density <- c(nnd_density,0)
                subject <- c(subject,x)
                
                par(mar = c(1, 1, 1, 1))
                plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),cex.axis = 0.8, cex.main = 0.8,main = x)
                
              } else if((i==j) == TRUE){
                g <- envelope(p1, Gest, nsim=sim, fix.n=TRUE, fix.marks=TRUE, savefuns=TRUE, nrank = m,verbose = FALSE)
                plot(g,legend=FALSE,xlab="", ylab="",cex.axis = 0.8, cex.main = 0.8,main = x)
                
                # Histogram of NND
                nn <- nndist(p1)
                average_nnd <- c(average_nnd,mean(nn))
                nnd_density <- c(nnd_density,nn)
                subject <- c(subject,rep(x,length(nn)))
                
                d <- distance
                if (max(g$r) < d){
                  dvalue <- max(g$r)
                } else {
                  dvalue <- distance 
                }
                index_value <- which.min(abs(g$r - dvalue))
                calculate_AUC(index = index_value, d = dvalue)
                
              } else {
                pm = ppp(c(p1x, p2x),c(p1y,p2y),window= window,marks=factor(c(rep(i,p1$n),rep(j,p2$n))))  
                g <- envelope(pm,Gcross, i=i, j=j, nsim = sim, fix.n=TRUE, fix.marks=TRUE, savefuns=TRUE, nrank = m, verbose = FALSE) # Gcross
                plot(g,legend=FALSE,xlab="", ylab="",cex.axis = 0.8, cex.main = 0.8,main = x)
                
                ## distance to nn
                nn <- nncross(p1,p2) # nncross
                nd2 <- as.data.frame(nn)
                average_nnd <- c(average_nnd,mean(nd2$dist))
                nnd_density <- c(nnd_density,nd2$dist)
                subject <- c(subject,rep(x,nrow(nd2)))
                
                d <- distance
                if (max(g$r) < d){
                  dvalue <- max(g$r)
                } else {
                  dvalue <- distance 
                }
                index_value <- which.min(abs(g$r - dvalue))
                calculate_AUC(index = index_value, d = dvalue)
              }
            }
            dev.off()
            
            ### plot histogram
            nd2 <- as.data.frame(list(nnd_density,subject),col.names = c("nn","sample_ID"))
            plot.new()
            pdf(paste0(diet,"_H-%03d.pdf"),width=11.4, height=6.5,onefile = FALSE)
            layout(matrix(c(1:28), ncol = 7 , byrow=TRUE))
            par(mfrow = c(4, 7),mgp=c(1,0.3,0),mar = c(2,1.5,1.5,1.5),oma = c(1, 1, 0, 0))
            pl <- unique(nd2$sample_ID)[2]
            for (pl in unique(nd2$sample_ID)){
              hd <- nd2 %>% filter(sample_ID == pl)
              hist(hd$nn, prob = TRUE, xlab = "distance",main = pl,cex.main=0.7,cex.lab=0.6,cex.axis=0.6,xlim=c(0,200))
              if (nrow(hd) > 1) {
                lines(density(hd$nn), xlim = c(0, 200))
              } else {}
            }
            dev.off()
          } 
          
          ### boxplot
          df_summary <- as.data.frame(list(sample_ID,p1n,p2n,AUC_cluster,dist_c,average_nnd,group_ID),col.names= c("ID",i,j,"AUC","dist_list","average_NND","group"))
          write.csv(df_summary, "df_summary.csv", row.names = FALSE)
          
          df_summary$AUC <- as.numeric(df_summary$AUC)
          df_summary <- na.omit(df_summary)
          
          stat.test2 <- suppressMessages(df_summary %>% na.omit() %>% 
                                           wilcox_test(AUC ~ group) %>%
                                           add_significance("p") %>% 
                                           add_xy_position())
          
          plot.new()
          filename <- "boxplot.pdf"
          ggplot(df_summary, aes(x=group,y=AUC)) + geom_boxplot(aes(fill = group),outlier.shape = NA, width = 0.4) + 
            geom_point(position=position_dodge(width=0.75),aes(group = group)) + theme_bw() + theme_classic() +
            stat_pvalue_manual(bracket.nudge.y = max(df_summary$AUC), stat.test2, label = "p.signif", tip.length = 0, hide.ns = TRUE, size = 10) +
            theme(legend.position = "none", axis.text.y = element_text(colour ="black",size = 13), axis.text.x = element_text(colour ="black",size = 13), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
            ylab("g area value")
          ggsave(filename, device = "pdf", width = 5, height = 7)
          dev.off()
        }
        
        
      } else {print("unique group not equal to 2")}
      
      print("DONE")
      
      
    } else {
      print("coordinate columns not integer")
    }
  })
}

shinyApp(ui = ui, server = server)