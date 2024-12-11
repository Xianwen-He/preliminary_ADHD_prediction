### This script aims to plot some plots 
### concerning the relationship between the DSM scores and the traits in interest.
### @Xianwen He, Nov 25, 2024

library(ggplot2)

# customized theme
mytheme <- theme(
  panel.background = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  axis.title = element_text(size = 14, hjust=0.5, family='serif'),
  axis.text = element_text(size = 12, family='serif'),  
  plot.title = element_text( # Center the title
    hjust = 0.5, size = 16, face = "bold", family='serif')
)
myblue <- 'cornflowerblue'

# read in the data
load('./Data/demographic.data.RData')

# ADH vs gender
adh.gender.plot <- ggplot(demographic.data%>%filter(!is.na(gender_factor)) , 
                          aes(x = gender_factor, y = DSM_Adh_Raw)) +
  geom_boxplot(fill = myblue, color = "black", alpha=0.8) +
  labs( 
    title = "ADH Scores by Gender",
    x = "Gender", y = "ADH Score") + 
  mytheme
ggsave(filename='./Data/img/adh_gender_plot.png', plot=adh.gender.plot,
       width=8, height = 6, dpi=200)

# Illicit substance
adh.illicit.plot <- ggplot(demographic.data%>%filter(!is.na(SSAGA_Times_Used_Illicits_binary)),
                          aes(x = SSAGA_Times_Used_Illicits_binary, y = DSM_Adh_Raw)) +
  geom_boxplot(fill = myblue, color = "black") +
  labs( 
    title = "ADH Scores by Illicit Substance Use",
    x = "Illicit Substance Use", y = "ADH Score") + 
  mytheme
ggsave(filename='./Data/img/adh_illicit_plot.png', plot=adh.illicit.plot,
       width=8, height = 6, dpi=200)

# Stimulant
adh.stimulant.plot <- ggplot(demographic.data%>%filter(!is.na(SSAGA_Times_Used_Stimulants_binary)),
                           aes(x = SSAGA_Times_Used_Stimulants_binary, y = DSM_Adh_Raw)) +
  geom_boxplot(fill = myblue, color = "black") +
  labs( 
    title = "ADH Scores by Stimulant Use",
    x = "Stimulant Use", y = "ADH Score") + 
  mytheme
ggsave(filename='./Data/img/adh_stimulant_plot.png', plot=adh.stimulant.plot,
       width=8, height = 6, dpi=200)

# Drunk
adh.alcohol.plot <- ggplot(demographic.data%>%filter(!is.na(SSAGA_Alc_12_Frq_Drk_binary)),
                             aes(x = SSAGA_Alc_12_Frq_Drk_binary, y = DSM_Adh_Raw)) +
  geom_boxplot(fill = myblue, color = "black") +
  labs( 
    title = "ADH Scores by Alcohol Use",
    x = "Whether Drunk", y = "ADH Score") + 
  mytheme
ggsave(filename='./Data/img/adh_alcohol_plot.png', plot=adh.alcohol.plot,
       width=8, height = 6, dpi=200)

# Conduct problems
adh.conduct.plot <- ggplot(demographic.data%>%filter(!is.na(SSAGA_ChildhoodConduct_binary)),
                           aes(x = SSAGA_ChildhoodConduct_binary, y = DSM_Adh_Raw)) +
  geom_boxplot(fill = myblue, color = "black") +
  labs( 
    title = "ADH Scores by Childhood Conduct Problem",
    x = "Problem Report", y = "ADH Score") + 
  mytheme
ggsave(filename='./Data/img/adh_conduct_plot.png', plot=adh.conduct.plot,
       width=8, height = 6, dpi=200)

# BMI
adh.bmi.plot <- ggplot(demographic.data, aes(x = BMI, y = DSM_Adh_Raw)) +
  geom_point(size = 2, color=myblue) +
  labs(
    title = "ADH Scores v.s. BMI",
    x = "BMI", y = "ADH Score") +
  mytheme
ggsave(filename='./Data/img/adh_bmi_plot.png', plot=adh.bmi.plot,
       width=8, height = 6, dpi=200)
