####### Multilevel Linguistic Analysis Across Modalities: Data Analysis #######

# In this script, we'll do our statistical analyses.

#############################################################################

##### 0. Preliminaries #####

# preliminaries
rm(list=ls())
library(lme4)
library(tidyverse)
library(lmerTest)
library(cowplot)
library(webshot)
setwd('~/GitHub/multilevel_linguistic_alignment_across_modalities/')

# thanks to https://rdrr.io/github/kylebutts/kfbmisc/src/R/ggplot_helpers.R
png_to_grob <- function(source) {
  
  temp <- png::readPNG(source)
  
  cli::cli_alert_info("Image has width = {dim(temp)[2]}px and height = {dim(temp)[1]}px")
  
  temp <- grid::rasterGrob(temp, interpolate = TRUE)
  return(temp)
}

##### 1. Data preparation #####

# load and prepare data
turn_df = read.csv('./data/turn_real-mlam.csv') %>%
  
  # extract the video from oen variable
  mutate(condition_info = gsub("_Transcript",
                               "",
                               condition_info,
                               ignore.case = TRUE)) %>%
  mutate(condition_info = gsub(".csv",
                               "",
                               condition_info,
                               ignore.case = TRUE)) %>%
  separate(condition_info,
           sep = "_",
           into = c("dyad",
                    "condition",
                    "conv_num",
                    "conv_type"),
           fill = "right") %>%
  mutate(conv_type = stringr::str_to_lower(conv_type)) %>%
  select(-X) %>%
  
  # fix some variable name location errors
  mutate(conv_type = ifelse(conv_num == "Arg",
                            "arg",
                            conv_type))

# let's plot what we have and see if we need to alter anything
hist(turn_df$lexical_lem1)
hist(log(turn_df$lexical_lem1))
hist(turn_df$lexical_lem2)
hist(log(turn_df$lexical_lem2))
hist(turn_df$syntax_penn_tok1)
hist(turn_df$syntax_penn_tok2)
hist(turn_df$cosine_semanticL)

# prepare for analysis
turn_df = turn_df %>% ungroup() %>%
  
  # log-transform where needed
  mutate(log_lexical_lem1 = log(turn_df$lexical_lem1),
         log_lexical_lem2 = log(turn_df$lexical_lem2)) %>%
  
  # transform INF to NA
  mutate(log_lexical_lem1 = ifelse(is.infinite(log_lexical_lem1),
                                   NA,
                                   log_lexical_lem1),
         log_lexical_lem2 = ifelse(is.infinite(log_lexical_lem2),
                                   NA,
                                   log_lexical_lem2)) %>%
  
  # recode for analyses
  mutate(conv_type = factor(conv_type, 
                            labels = c("Aff.", 
                                       "Arg.", 
                                       "Coop."))) %>%
  mutate(condition = factor(condition, 
                            labels = c("FF Lab", 
                                       "VC Remote", 
                                       "VC Lab"))) %>%
  mutate(condition = as.factor(condition),
         conv_type = as.factor(conv_type),
         dyad = as.factor(dyad))
contrasts(turn_df$condition) = contr.treatment(3)
contrasts(turn_df$conv_type) = contr.treatment(3)

##### 2. Unigram lexical alignment #####

# analysis: unigram lexical alignment by time, condition, and conversation type
model_lexical_unigram_alignment = lmer(log_lexical_lem1 ~ time + condition * conv_type +
                                         (1 | dyad),
                                       data = turn_df)
summary(model_lexical_unigram_alignment)

# table: unigram lexical alignment results
table_lexical_unigram_alignment = sjPlot::tab_model(model_lexical_unigram_alignment,
                                                    dv.labels = c("Lexical Unigram Alignment"),
                                                    pred.labels = c("Intercept",
                                                                    "Time",
                                                                    "VC Remote",
                                                                    "VC Lab",
                                                                    "Arg.",
                                                                    "Coop.",
                                                                    "VC Remote x Arg.",
                                                                    "VC Lab x Arg.",
                                                                    "VC Remote x Coop.",
                                                                    "VC Lab x Arg."),
                                                    file = "./tables/table_lexical_unigram_alignment.html")
table_lexical_unigram_alignment

# plot: unigram lexical alignment
plot_unigram_lexical_alignment = ggplot(turn_df,
                                        aes(y = log_lexical_lem1,
                                            x = conv_type,
                                            color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  theme(legend.position = "none") +
  geom_jitter(width=.2,
              alpha=.1) + 
  xlab("Conversation type")  +
  ylab("Log Unigram Lexical Alignment") +
  ggtitle("Lexical unigram alignment by\ncondition and conversation type")
ggsave(filename = paste0('./figures/plot_unigram_lexical_alignment-mlam.png'),
       plot = plot_unigram_lexical_alignment,
       height = 4,
       width = 3.05,
       units = "in")

##### 3. Bigram lexical alignment #####

# analysis: bigram lexical alignment by time, condition, and conversation type
model_lexical_bigram_alignment = lmer(log_lexical_lem2 ~ time + condition * conv_type +
                                        (1 + conv_type | dyad),
                                      data = turn_df)
summary(model_lexical_bigram_alignment)

# table: bigram lexical alignment results
table_lexical_bigram_alignment = sjPlot::tab_model(model_lexical_bigram_alignment,
                                                   dv.labels = c("Lexical Bigram Alignment"),
                                                   pred.labels = c("Intercept",
                                                                   "Time",
                                                                   "VC Remote",
                                                                   "VC Lab",
                                                                   "Arg.",
                                                                   "Coop.",
                                                                   "VC Remote x Arg.",
                                                                   "VC Lab x Arg.",
                                                                   "VC Remote x Coop.",
                                                                   "VC Lab x Arg."),
                                                   file = "./tables/table_lexical_bigram_alignment.html")
table_lexical_bigram_alignment

##### 3. Unigram syntactic alignment #####

# analysis: unigram syntactic alignment by time, condition, and conversation type
model_syntactic_unigram_alignment = lmer(syntax_penn_lem1 ~ time + condition * conv_type +
                                           (1 + conv_type | dyad),
                                         data = turn_df)
summary(model_syntactic_unigram_alignment)

# table: unigram syntactic alignment results
table_syntactic_unigram_alignment = sjPlot::tab_model(model_syntactic_unigram_alignment,
                                                      dv.labels = c("Syntactic Unigram Alignment"),
                                                      pred.labels = c("Intercept",
                                                                      "Time",
                                                                      "VC Remote",
                                                                      "VC Lab",
                                                                      "Arg.",
                                                                      "Coop.",
                                                                      "VC Remote x Arg.",
                                                                      "VC Lab x Arg.",
                                                                      "VC Remote x Coop.",
                                                                      "VC Lab x Arg."),
                                                      file = "./tables/table_syntactic_unigram_alignment.html")
table_syntactic_unigram_alignment

# plot: unigram syntactic alignment
plot_unigram_syntactic_alignment = ggplot(turn_df,
                                          aes(y = syntax_penn_lem1,
                                              x = conv_type,
                                              color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  theme(legend.position = "bottom") +
  geom_jitter(width=.2,
              alpha=.3) + 
  xlab("Conversation type")  +
  ylab("Unigram Syntactic Alignment") +
  ggtitle("Syntactic unigram alignment\nby condition and conversation type")
ggsave(filename = paste0('./figures/plot_unigram_syntactic_alignment-mlam.png'),
       plot = plot_unigram_syntactic_alignment,
       height = 6,
       width = 4,
       units = "in")

##### 4. Bigram syntactic alignment #####

# analysis: bigram syntactic alignment by time, condition, and conversation type
model_syntactic_bigram_alignment = lmer(syntax_penn_lem2 ~ time + condition * conv_type +
                                          (1 + conv_type | dyad),
                                        data = turn_df)
summary(model_syntactic_bigram_alignment)

# table: bigram syntactic alignment results
table_syntactic_bigram_alignment = sjPlot::tab_model(model_syntactic_bigram_alignment,
                                                     dv.labels = c("Syntactic Bigram Alignment"),
                                                     pred.labels = c("Intercept",
                                                                     "Time",
                                                                     "VC Remote",
                                                                     "VC Lab",
                                                                     "Arg.",
                                                                     "Coop.",
                                                                     "VC Remote x Arg.",
                                                                     "VC Lab x Arg.",
                                                                     "VC Remote x Coop.",
                                                                     "VC Lab x Arg."),
                                                     file = "./tables/table_syntactic_bigram_alignment.html")
table_syntactic_bigram_alignment

# plot: bigram syntactic alignment
plot_bigram_syntactic_alignment = ggplot(turn_df,
                                         aes(y = syntax_penn_lem2,
                                             x = conv_type,
                                             color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  theme(legend.position = "bottom") +
  geom_jitter(width=.2,
              alpha=.3) + 
  xlab("Conversation type")  +
  ylab("Bigram Syntactic Alignment") +
  ggtitle("Syntactic bigram alignment\nby condition and conversation type")
ggsave(filename = paste0('./figures/plot_bigram_syntactic_alignment-mlam.png'),
       plot = plot_bigram_syntactic_alignment,
       height = 6,
       width = 4,
       units = "in")

##### 5. Conceptual alignment #####

# analysis: conceptual alignment by time, condition, and conversation type
model_conceptual_alignment = lmer(cosine_semanticL ~ time + condition * conv_type +
                                    (1 | dyad),
                                  data = turn_df)
summary(model_conceptual_alignment)

# table: conceptual alignment results
table_conceptual_alignment = sjPlot::tab_model(model_conceptual_alignment,
                                               dv.labels = c("Conceptual Alignment"),
                                               pred.labels = c("Intercept",
                                                               "Time",
                                                               "VC Remote",
                                                               "VC Lab",
                                                               "Arg.",
                                                               "Coop.",
                                                               "VC Remote x Arg.",
                                                               "VC Lab x Arg.",
                                                               "VC Remote x Coop.",
                                                               "VC Lab x Arg."),
                                               file = "./tables/table_conceptual_alignment.html")
table_conceptual_alignment

# plot: conceptual alignment
plot_conceptual_alignment = ggplot(turn_df,
                                   aes(y = cosine_semanticL,
                                       x = conv_type,
                                       color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  theme(legend.position = "bottom") +
  geom_jitter(width=.2,
              alpha=.3) + 
  xlab("Conversation type")  +
  ylab("Conceptual Alignment") +
  ggtitle("Conceptual alignment\nby condition and conversation type")
ggsave(filename = paste0('./figures/plot_conceptual_alignment-mlam.png'),
       plot = plot_conceptual_alignment,
       height = 6,
       width = 4,
       units = "in")

##### 6. Multi-paneled figure for CogSci 2023 submission #####

# plot: unigram lexical alignment
plot_subpanel_unigram_lexical_alignment = ggplot(turn_df,
                                                 aes(y = log_lexical_lem1,
                                                     x = conv_type,
                                                     color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  theme(plot.margin = margin(6, 0, 6, 0))+ 
  geom_jitter(width=.2,
              alpha=.1) + 
  xlab("")  +
  ylab("Log Alignment") +
  ggtitle("Lexical unigram")

# plot: unigram syntactic alignment
plot_subpanel_unigram_syntactic_alignment = ggplot(turn_df,
                                                   aes(y = syntax_penn_lem1,
                                                       x = conv_type,
                                                       color = conv_type)) +
  geom_violin() +
  facet_grid(cols = ,
             rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conversation Type")+
  geom_jitter(width=.2,
              alpha=.1) + 
  theme(plot.margin = margin(6, 0, 6, 0))+ 
  coord_cartesian(ylim = c(0,1))+
  xlab("")  +
  ylab("Syntactic alignment") +
  ggtitle("Unigram")

# plot: unigram and bigram syntactic alignment
plot_subpanel_bigram_syntactic_alignment = ggplot(turn_df,
                                                  aes(y = syntax_penn_lem2,
                                                      x = conv_type,
                                                      color = conv_type)) +
  geom_violin() +
  facet_grid(rows = vars(condition)) +
  scale_color_manual(values = c("blue", "red","darkgray"),
                     labels = c("Aff.",
                                "Arg.",
                                "Coop."),
                     name = "Conv. Type")+
  theme(plot.margin = margin(6, 0, 6, 0))+ 
  coord_cartesian(ylim = c(0,1))+
  geom_jitter(width=.2,
              alpha=.1) + 
  xlab("")  +
  ylab("") +
  ggtitle("Bigram")

# create subpanels and legend
panel_unigram_syntactic_alignment = cowplot::plot_grid(
  plot_subpanel_unigram_syntactic_alignment + theme(legend.position = "none"),
  labels = c("A")
)
panel_bigram_syntactic_alignment = cowplot::plot_grid(
  plot_subpanel_bigram_syntactic_alignment + theme(legend.position = "none"),
  labels = c("B")
)
plot_alignment_legend = get_legend(plot_subpanel_bigram_syntactic_alignment + 
                                     theme(legend.position = "bottom", 
                                           legend.box.margin = margin(0, 0, 0, 0)))
plot_alignment_data = cowplot::plot_grid(panel_unigram_syntactic_alignment, 
                                         panel_bigram_syntactic_alignment, 
                                         nrow=1)
# create and add title
plot_alignment_title = ggdraw() + 
  draw_label("Syntactic alignment by\ncondition and conversation type", 
             fontface = 'bold', x = 0, hjust = 0) 
cogsci_fig_2 = cowplot::plot_grid(plot_alignment_title,
                                  plot_alignment_data,
                                  # rel_heights = c(0.1, 1, .1),
                                  rel_heights = c(0.2, 1),
                                  ncol=1)
cogsci_fig_2

# save to file
ggsave(filename = paste0('./figures/cogsci_fig_2-mlam.png'),
       plot = cogsci_fig_2,
       height = 4,
       width = 4,
       units = "in")
 
##### 7. Multi-paneled table for CogSci 2023 submission #####

# subtable: lexical unigram alignment
table_subtable_lexical_unigram_alignment = sjPlot::tab_model(model_lexical_unigram_alignment,
                                                             dv.labels = c("Lexical Unigram Alignment"),
                                                             pred.labels = c("Intercept",
                                                                             "Time",
                                                                             "VC Remote",
                                                                             "VC Lab",
                                                                             "Arg.",
                                                                             "Coop.",
                                                                             "VC Remote x Arg.",
                                                                             "VC Lab x Arg.",
                                                                             "VC Remote x Coop.",
                                                                             "VC Lab x Arg."),
                                                             show.re.var = FALSE,
                                                             show.obs = FALSE,
                                                             show.icc = FALSE,
                                                             show.ngroups = FALSE,
                                                             file = "./tables/table_subtable_lexical_unigram_alignment.html")
table_subtable_lexical_unigram_alignment
webshot("./tables/table_subtable_lexical_unigram_alignment.html",
        "./tables/table_subtable_lexical_unigram_alignment.png",
        selector = 'table') %>%
  resize("150%")
grob_subtable_lexical_unigram_alignment = png_to_grob("./tables/table_subtable_lexical_unigram_alignment.png")

# subtable: lexical bigram alignment
table_subtable_lexical_bigram_alignment = sjPlot::tab_model(model_lexical_bigram_alignment,
                                                            dv.labels = c("Lexical Bigram Alignment"),
                                                            pred.labels = c("Intercept",
                                                                            "Time",
                                                                            "VC Remote",
                                                                            "VC Lab",
                                                                            "Arg.",
                                                                            "Coop.",
                                                                            "VC Remote x Arg.",
                                                                            "VC Lab x Arg.",
                                                                            "VC Remote x Coop.",
                                                                            "VC Lab x Arg."),
                                                            show.re.var = FALSE,
                                                            show.icc = FALSE,
                                                            show.obs = FALSE,
                                                            show.ngroups = FALSE,
                                                            file = "./tables/table_subtable_lexical_bigram_alignment.html")
table_subtable_lexical_bigram_alignment
webshot("./tables/table_subtable_lexical_bigram_alignment.html",
        "./tables/table_subtable_lexical_bigram_alignment.png",
        selector = 'table') %>%
  resize("150%")
grob_subtable_lexical_bigram_alignment = png_to_grob("./tables/table_subtable_lexical_bigram_alignment.png")

# bind grobs in cowplot
table_space = grid::nullGrob()
table_lexical_alignment_data = cowplot::plot_grid(grob_subtable_lexical_unigram_alignment,
                                                  table_space,
                                                  grob_subtable_lexical_bigram_alignment,
                                                  rel_widths = c(1,.01,1),
                                                  nrow=1)

# create and add title
table_alignment_title = ggdraw() + 
  draw_label("Models predicting lexical alignment", 
             fontface = 'bold', x = 0, hjust = 0) 
cogsci_table_1 = cowplot::plot_grid(table_alignment_title,
                                    table_lexical_alignment_data,
                                    rel_heights = c(0.1, 1),
                                    ncol=1)
cogsci_table_1

# save to file
ggsave(filename = paste0('./tables/cogsci_table_1-mlam.png'),
       plot = cogsci_table_1,
       height = 3,
       width = 6,
       units = "in")
 