require(flexplot)
require(tidyverse)
require(lme4)
require(party)
library(dplyr)
mentalhealth= read.csv("data/W1mentalhealth.csv")
head(mentalhealth)
personality=read.csv("data/W1personality .csv")
head(personality)
physicalsensory=read.csv("data/W1physicalsensory.csv")
socialinteraction=read.csv("data/W1socialinteraction.csv")
wellbeing=read.csv("data/W1wellbeing.csv")
?left_join
#left join keeps all of x and adds onto it y 
view(mentalhealth)
#137 columns 
view(personality)
#105 columns 
combined1=left_join(mentalhealth,personality)
view(combined1)
#236 columns 
combined2=left_join(combined1, physicalsensory)
view(combined2)
#297 columns 
combined3= left_join(combined2, socialinteraction)
view(combined3)
#394 columns 
d_original=combined4=left_join(combined3,wellbeing)
view(d_original)
#440 columns 
head(d_original)
d= d_original %>% 
  mutate(sum_cati=rowSums(select(., contains("cati"))))
head(d)
flexplot(sum_cati~1, data=d)
d$sum_cati
#Something not working


#grep function is like select but prints the column names 
#Try a different measure 







