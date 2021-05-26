########################################################################
# Tess Linden
# May 2021 
#
# make_prot_fig.R
#
# defines the function make_prot_fig(), which takes as input: (1) a data 
# frame of proteins and their domains with their respective start and 
# stop locations; and (2) a data frame that assigns each domain a color,
# and generates a diagram of the proteins and their domains. 
#
# when called with the argument legend = TRUE, the function will plot 
# the diagram's legend instead of the diagram itself.
#
# the input df "colors" should have at least the following columns:
# domain, color
# (colors can be in the form of hex codes)
#
# the input df "data" should have at least the following columns:
# protein_id, domain, start, stop
# e.g. the row protein_x, Kinase, 5, 55 means protein_x has a Kinase 
# domain from amino acid #5 to amino acid #55.
#
# also, for each protein in the "data" df, there should exist a row 
# that defines the length of that protein, formatted as follows:
# [protein id], protein, 1, [protein length].
#
#
#
########################################################################




make_prot_fig <- function(data, colors, scalebar=500, title=NULL, legend = FALSE)
{
  
  # TODO: Throw error if any domains overlap

  pr_height <- 1
  dom_height <- 10
  sep <- 10
  
  # Print warning if any protein-length rows are missing
  for(protein_id in unique(data$protein_id))
  {
    sub <- data[data$protein_id == protein_id,]
    if(!('protein' %in% sub$domain))
    {
      print(paste0('WARNING: Protein \'',protein_id,'\' is missing a protein-length row.'))
    }
  }
  
  # Print warning if any domains are missing color assignments
  if(sum(!(unique(data$domain[data$domain != 'protein']) %in% colors$domain)) != 0)
  {
    print('WARNING: One or more domains is missing a color assignment.')
  }
  
  # make sure protein ids and domains are encoded as strings, not factors
  data$protein_id <- as.character(data$protein_id)
  data$domain <- as.character(data$domain)
  colors$domain <- as.character(colors$domain)
  colors$color <- as.character(colors$color)
  
  # get rid of any commas in numbers 
  data$start <- as.numeric(gsub(",", "", data$start))
  data$stop <- as.numeric(gsub(",", "", data$stop))
  
  # add length column
  data$length <- data$stop - data$start + 1
  
  # if legend = TRUE, plot legend
  if(legend)
  {
    par(mar = c(5.1, 4.1, 4.1, 2.1)) # resets margins to defaults
    plot(x=0, xlab='', ylab='', type='n', main = title, xaxt='n', yaxt='n', frame.plot=FALSE)
    legend <- colors[colors$domain %in% data$domain,]
    legend(x=0, legend=legend$domain, fill=legend$color, cex=0.8)
  }
  # else, plot proteins
  else
  {
    # set plot's left margin width based on length of longest protein id
    par(mai=c(1, 0.5+max(strwidth(data$protein_id,units='inches')), 1, 1), xpd=TRUE) #syntax: c(bottom, left, top, right)
    
    # set x and y ranges
    num_proteins <- length(unique(data$protein_id))
    ylim <- c(0, num_proteins*(sep+dom_height)+sep*2.5)
    xlim <- c(-100, max(data$length)+100)
    if(xlim[2] <= scalebar) { xlim[2] <- scalebar+1 }
    
    # make blank plot
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', xaxt='n',yaxt='n', frame.plot=FALSE)
    text(x = 0, y = ylim[2], adj = 0, labels = title, font = 2, cex = 1.25)
    
    # plot scalebar
    axis(side = 1, at = c(0,scalebar), labels = FALSE)
    text(x = scalebar/2, y = 0, adj = 0.5, labels=paste(scalebar,'aa'))
    
    
    # iterate over proteins in data set
    y <- max(ylim)-sep/2
    for(protein in unique(data$protein_id))
    {
      curr_protein <- data[data$protein_id == protein,]
      # calculate the y coordinate of this protein 
      y <- y - sep - dom_height
      # plot the protein backbone
      rect(1, y-pr_height/2, max(curr_protein$length), y+pr_height/2, col='#000000')#, lwd = 1)
      
      # add protein label
      text(x=-100,y=y,labels=protein, adj=1) #adj=0 means left-aligned, 1 means right-aligned
      
      # plot domains onto backbone
      domains <- unique(curr_protein$domain)[unique(curr_protein$domain) != 'protein']
      for (dom in domains)
      {
        curr_domain <- curr_protein[curr_protein$domain == dom,]
        for (i in 1:length(curr_domain[,1]))
        {
          rect(curr_domain$start[i], y-dom_height/2, curr_domain$stop[i], y+dom_height/2, 
               col = colors$color[colors$domain == dom]) #col = colors[j+1])#, lwd = 1)
        }
      }
    }
  }
} 
# end definition of make_prot_fig()




#
# sample usage 
# 
setwd("/Users/tesslinden/Dropbox/King Lab/KL Notebook/docker-singularity-etc/scripts/for\ sharing")
my_data = read.csv('prot_fig_sample_data.csv', na.strings = '', stringsAsFactors = FALSE)
my_colors <- data.frame(domain = unique(my_data$domain[my_data$domain != 'protein']),
                        color = c('#56B4E9','#E69F00','#009E73','#F0E442'))
# plot diagrams
make_prot_fig(my_data, my_colors, title = 'sample proteins')
# plot legend
make_prot_fig(my_data, my_colors, title = 'sample proteins', legend = TRUE)
