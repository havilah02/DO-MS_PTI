init <- function() {
  
  type <- 'plot'
  box_title <- 'Number of Confident Identifications'
  help_text <- 'Plotting the number of peptides identified at each given confidence level.'
  source_file <- 'evidence'
  
  .validate <- function(data, input) {
    validate(need(data()[['evidence']], paste0('Upload evidence.txt')))
  }
  
  .plotdata <- function(data, input) {
    plotdata <- data()[['evidence']][,c('Raw.file', 'PEP','Type')]
    plotdata <- plotdata %>% dplyr::filter(Type != "MULTI-MATCH")
    plotdata <- plotdata %>% dplyr::select('Raw.file', 'PEP')
    
    # creating a sequence of logarithmically spaced bins for the PEP values
    # build log10 PEP vector
    min_pep <- max(c(min(plotdata$PEP, na.rm = TRUE), 1e-5)) # making sure the minimum is 1e-5
    max_pep <- max(plotdata$PEP, na.rm = TRUE)
    peps <- seq(log10(min_pep), log10(max_pep), length.out = 500) # create 500 bins
    peps <- c(log10(.Machine$double.xmin), peps) # add machine minimum
    # peps <- seq(log10(max(c(min(plotdata$PEP)), 1e-5)), log10(max(plotdata$PEP)), length.out=500)
    # peps <- c(log10(.Machine$double.xmin), peps)
    
    plotdata <- plotdata %>%  #error
      dplyr::mutate(bin=findInterval(PEP, 10**peps, rightmost.closed = TRUE)) %>%
      dplyr::group_by(Raw.file, bin) %>%
      dplyr::summarise(n=dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(cy=cumsum(n), # cumulative number of peptides
                    pep = ifelse(bin <= length(peps), 10**peps[bin], NA) #bin to pep conversion
      ) %>%
      dplyr::filter(!is.na(pep)) # remove any rows with invalid bin references
    #pep=10**peps[bin])
    return(plotdata)
  }
  
  .plot <- function(data, input) {
    .validate(data, input)
    plotdata <- .plotdata(data, input)  #error
    
    validate(need((nrow(plotdata) > 1), paste0('No Rows selected')))
    
    # rank the experiments by cumulative ID
    rank_info <- plotdata %>%
      dplyr::group_by(Raw.file) %>%
      dplyr::summarise(max_cy = max(cy, na.rm = TRUE)) %>%
      dplyr::arrange(desc(max_cy)) %>%
      dplyr::mutate(rank_ord = row_number())
    
    plotdata <- plotdata %>%
      dplyr::left_join(rank_info, by = "Raw.file")
    
    cc <- scales::seq_gradient_pal('red', 'blue', 'Lab')(seq(0, 1, length.out = nrow(rank_info)))
    
    ggplot(plotdata, aes(x=pep, color=factor(rank_ord), y=cy, group=Raw.file)) + 
      geom_line(size = input$figure_line_width) +
      scale_colour_manual(name='Experiment', values=cc, labels=rank_info$Raw.file) +
      coord_flip() + 
      scale_x_log10(limits=c(.00009,.1), breaks=c(.0001,.001,.01,.1), 
                    labels=scales::trans_format('log10', scales::math_format(10^.x))) + 
      theme_base(input=input) +
      theme(legend.position='right') + 
      theme(legend.key=element_rect(fill='white')) +
      xlab('PEP') + ylab('Number of IDs') 
    
  }
  
  return(list(
    type=type,
    box_title=box_title,
    help_text=help_text,
    source_file=source_file,
    validate_func=.validate,
    plotdata_func=.plotdata,
    plot_func=.plot,
    box_width=12, # bootstrap column units
    plot_height=500, # pixels
    report_plot_width=7, # inches
    report_plot_height=5 # inches
  ))
}
