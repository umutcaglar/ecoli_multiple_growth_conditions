# This is the source of the bug 
# We are looking at 
# a) cophenetic_distance_RNA_main.R lines 191-200
# b) cophenetic_distance_protein_main.R  190-199

# correct way of doing it
main_growthPhaseF %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary_correct_1

main_growthPhaseF %>%
  dplyr::group_by(iteration) %>%
  dplyr::mutate(overall_Mean=unique(overal_Mean))%>%
  dplyr::group_by() %>%
  dplyr::summarize(overall_fakeMean=mean(overall_Mean),
                   overall_fakeStd=sd(overall_Mean))->FakeSummary_correct_2

FakeSummary_correct<-cbind(FakeSummary_correct_1,FakeSummary_correct_2)



# the way I did it
main_growthPhaseF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal), # problematic lines 1
                overall_fakeStd=sd(meanVal))%>% # problematic lines 2
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary_false


# the alternative way
main_growthPhaseF %>%
  dplyr::mutate(overall_fakeMean=mean(overal_Mean),
                overall_fakeStd=sd(overal_Mean))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary_alternative