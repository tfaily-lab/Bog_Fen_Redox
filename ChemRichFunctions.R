chemRICHFix <- function(x){
  #Data Check:
  if(sum(colnames(x) %in% c('compound_name', 'effect_size', 'pvalue', 'set')) != 4){
    stop('Data Incorrect: Make sure column names are compound_name, effect_size, pvalue, and set')
  }
  
  x_out <- x %>%
    #Decide the direction that metabolites are moving 
    mutate(direction = ifelse(effect_size > 0, 'up', 'down')) %>%
    #Non-significant metabolites are not moving
    mutate(direction = ifelse(pvalue > 0.05, 'no_change', direction)) %>%
    #Conver effect size for downward moving metabolites into 1/efs
    mutate(efs = ifelse(effect_size < 0, 1/abs(effect_size), abs(effect_size))) %>%
    #Non-moving metabolites have no effect size
    mutate(efs = ifelse(pvalue > 0.05, 1, efs)) %>%
    #Drop the NAs from the dataset
    filter(!is.na(set)) %>%
    group_by(set) %>%
    #Remove any compound classes which have no significant differences
    filter(!all(pvalue > 0.05)) %>%
    #Remove any classes with less than 2 compounds
    filter(n() > 2) %>%
    nest() %>%
    #Run the KS-Test essentially asking the question, what is the likelihood the pvalue distribution is uniform for this metabolite?
    mutate(ptest = map_dbl(data, function(x) suppressWarnings(ks.test(x$pvalue, 'punif', alternative = 'greater')$p.value))) %>%
    #Correct p-values rounded to 0
    modify_at('ptest', ~ifelse(.x == 0, 2.2e-20, .x)) %>%
    #unnest('data') 
    #Calculate ration of significant to all
    mutate(altrat = map_dbl(data, function(x) sum(x$pvalue < 0.05)/length(x$pvalue))) %>%
    #Calculate the ratio of metabolites which were increased
    mutate(uprat = map_dbl(data, function(x) sum(x$pvalue < 0.05 & x$direction == 'up')/sum(x$pvalue < 0.05))) %>% 
    #Calculate the number of metaboltes in the group
    mutate(size = map_dbl(data, function(x) nrow(x))) %>%
    #Calculate the number of significant metabolties
    mutate(num_altered = map_dbl(data, function(x) sum(x$pvalue < 0.05))) %>%
    #How many went up?
    mutate(num_up = map_dbl(data, function(x) sum(x$direction == 'up' ))) %>%
    #How many went down?
    mutate(num_down = map_dbl(data, function(x) sum(x$direction == 'down' ))) %>%
    #Ungroup to prevent weird things from happengin
    ungroup() %>%
    #Now adjust p-values
    modify_at('ptest', ~p.adjust(.x, method = 'fdr')) %>%
    #Create a random order for chemical classes
    mutate(order = sample(1:nrow(.), replace = F))
  return(x_out)
}

suppressMW <- function(x){
  suppressMessages(suppressWarnings(x))
}

plotChemRich <- function(CR_output, interactive = FALSE, colors = c('blue', 'red'), nset = 5){
  p <- suppressMW(ggplot(CR_output, aes(x = order, y = -log(ptest), size = size, color = uprat)) +
                    geom_point(aes(text = paste0('Class: ', set, '\n',
                                                 'K-S P-Val: ', ptest, '\n',
                                                 'Up Expressed: ', num_up, '\n',
                                                 'Down Expressed: ', num_down, '\n',
                                                 'Size: ', size))) +
                    scale_color_gradient(low = colors[1], high = colors[2], limits = c(0,1)) +
                    scale_size(range = c(5,30)) +
                    ggrepel::geom_label_repel(aes(label = set), color = 'gray20', data = subset(CR_output, size > nset), size = 4, force = 5) +
                    theme(legend.position = 'none'))
  if(interactive){
    suppressMW(plotly::ggplotly(p, tooltip = 'text'))
  }else{
    suppressMW(plot(p))
  }
} 

plotChemRich_custom <- function(CR_output, interactive = FALSE, colors = c('blue', 'white', 'red'), sets, 
                                text_size = 4, binary = FALSE, sig_line = 0.05){
  
  CR_output <- CR_output %>%
    mutate(Expression_Category = case_when(
      uprat <= 0.45 ~ "Enriched in Fen",
      uprat >= 0.65 ~ "Enriched in Bog",
      TRUE ~ "Enriched in Neither"
    ),
    Expression_Category = factor(Expression_Category, 
                                 levels = c("Enriched in Fen", "Enriched in Bog", "Enriched in Neither")))
  
  p <- suppressMW(
    ggplot(CR_output, aes(x = order, y = -log(ptest), size = size, color = Expression_Category)) +
      geom_hline(yintercept = -log(sig_line), linetype = 'dashed') +
      geom_point(aes(text = paste0('Class: ', set, '\n',
                                   'K-S P-Val: ', ptest, '\n',
                                   'Up Expressed: ', num_up, '\n',
                                   'Down Expressed: ', num_down, '\n',
                                   'Size: ', size))) +
      scale_color_manual(
        values = c("Enriched in Fen" = colors[1], "Enriched in Neither" = colors[2], "Enriched in Bog" = colors[3]),
        name = NULL  # Removes legend title
      ) +
      scale_size(range = c(5,30), guide = 'none') +
      ggrepel::geom_label_repel(aes(label = set), 
                                color = 'gray20', 
                                data = subset(CR_output, set %in% sets), 
                                size = text_size, force = 5) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(
        legend.position = c(0.85, 0.85), # Legend inside the top-right corner
        legend.background = element_rect(fill = alpha('white', 0.6), color = NA), # Optional: Transparent legend box
        legend.title = element_blank()
      )
  )  # <- this parenthesis closes ggplot()
  
  if(interactive){
    suppressMW(plotly::ggplotly(p, tooltip = 'text'))
  } else {
    suppressMW(plot(p))
  }
}



