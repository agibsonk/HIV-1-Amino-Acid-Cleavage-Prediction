# HIV-1 Amino Acid Cleavage Prediction 
Using SVC and kNN learning methods to predict whether or not HIV-1 polymerase will cleave any given 8 long amino acid chain.
Trained using a dataset of 1625 amino acid chains. Other datasets are available as well in sub-folder.
Analysis of dataset was also done to determine data distribution, such as checking what amino acid most chains primariy consisted of, as well as examining relationship between an amino acid being cleaved and whether the amino acid was composed of many or primarily one amino acid type.
The models were ~60% accurate, which is a marginal increase over a coin flip of whether the chain will be cleaved
Further improvements that could be made include: using cnn model to potentially improve accuracy, combining datasets to improve data diversity, and properly normalizing data.
