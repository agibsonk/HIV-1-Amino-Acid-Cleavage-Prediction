import csv
import numpy as np
from sklearn import svm
import scipy.stats
import matplotlib.pyplot as plt

##Classifiers:
##Cleaved = +1
##Not Cleaved = -1


def load_HIV_csv(filename):

    reader = csv.reader(open(filename, 'r'))
    data_list = []

    for row in reader:
        data_list.append(row)

    return data_list
	
def organize_data(data): ##Assumed data will come in as ['Letters','value'] and both are strings
	
	amino_acids = []
	observations = []
	
	
	for row in range(len(data)):
		
		amino_acids = [] ##Reseting the amino acid list to being empty
		for letter in range(len(data[row][0])):
			amino_acids += data[row][0][letter] ##Generating a list that is the letters of each amino acid
		observations += [amino_acids] ##Need to treat each amino acid list as its own list, not concatenate if there was no brackets
	
	
	for list in range(len(observations)): ##Want to convert the letters in list to their numerical counterparts
		for letter in range(len(observations[list])): 
			observations[list][letter] = ord(observations[list][letter])

	
	observations = np.array(observations)##Converting list of data to an array
	
	return(observations)
	
def organize_classifiers(data): ##Want to build a array of n X 1 of the classifiers from the data
	classifiers = []
	for row in range(len(data)):
		classifiers = classifiers + [data[row][1]]
	
	classifiers = [int(i) for i in classifiers] ##Converting the strings in classifiers to integers
	classifiers = np.array(classifiers) ##Converting it to an array
	return classifiers
	
def load_data(filename): ##Incorperates previous functions to load and process data from csv -> 2 np 
	"""
	Loads amino acid sequences and if they are cut from specified CSV file and stores it in 2 numpy arrays, one for the data (sequences) and another for the labels (if they are cut or not) NOTE: Amino acid sequences have been converted from letters to numerical values using python's ord() function
	:param filename: CSV file containing one row of amino acid sequences and another row of either +1 or -1 depending if cut or not
	:return: 2 numpy arrays, array 1 is sequence data, array 2 is labels
	"""

	working_data = load_HIV_csv(filename)##Importing CSV file
	clean_data = organize_data(working_data)##Making the data row of CSV list into array specifying each letter in amino acid sequence
	classifiers = organize_classifiers(working_data) ##Making the classifier row of CSV into array of all classifiers for the data
	
	return clean_data,classifiers


def viz_1(data,labels): ##Want to see most common letter (mode) that results in a cut
	"""
	Makes a bar graph that is determining the "mode" amino acid of each amino acid chain (ie.whichever amino acid is most common in each sequence) and plotting the total amount (regardless of if the sequence is cleaved or not) of these along with how many are cut overtop, showing a distribution of the range of data we are testing
	:param data: A numpy array that contains 8 rows of data, each row a numerical equivalent of the amino acid letter
	:param labels: A numpy array that is only 1 dimensional and specifies if the amino acid chain with the same index in the data parameter was cut or not (+1 = cut , -1= not cut)
	
	:return: A histogram that plots the number of chains that had a mode amino acid and overlays another bar that is the number of those modes that were cleaved
	"""


	modes = []
	frequency_of_cut_and_uncut = {} ##Making a dict showing how many letter modes there are
	frequency_of_cut = {}  ##Making a dict showing how many of these letter modes there are that get cut
	
	for i in range(len(data)):
		modes += list(scipy.stats.mode(data[i])[0]) ##Recording the mode of each observation (most common letter after/ord).
	
	for i in range(len(modes)):
		modes[i] = chr(modes[i]) ##Converts numbers to letters. 
	
	
	##Modes is now a list of the mode letter of each datapoint. Need to make a dictionary that increases in value.
	
	for unique in set(modes): ##Want to get the unique modes as 0 in both dictionaries
		frequency_of_cut[unique] = 0
		frequency_of_cut_and_uncut[unique] = 0
		
		
	for i in range(len(modes)): ##Making a dictionary that maps the mode letter of sequence -> How many times this resulted in a cut
		if labels[i] == 1:
			frequency_of_cut[modes[i]] += 1 ##Adding one only if it gets cut
		frequency_of_cut_and_uncut[modes[i]] +=1 ##Adding regardless if cut or uncut
	
	p1 = plt.bar(*zip(*frequency_of_cut_and_uncut.items()))
	p2 = plt.bar(*zip(*frequency_of_cut.items())) ##Plots a bar chart from dictionary
	plt.title('Amino Acid Chains Cleaved by HIV-1 Protease According to Majority Amino Acid')
	plt.ylabel('Occurences')
	plt.xlabel('Majority Amino Acid in Chain')
	plt.legend((p1[0],p2[0]),('Total Amount of Chains Tested','Portion of Chains Cleaved'))
	
	plt.show()


def viz_2(data,labels): ##Want a visualization of what percentage is cleaved depending on the portion of each amino acid that make up the chains
	"""
	Creates a scatterplot that plots each chain depending on what portion of the chain is made up by an amino acid, and what percentage of these chains were cut, and does this for all 20 amino acids. Since each chain is 8 amino acids long, there will be 20 points plotted for each portion. May show if having more of one amino acid in a chain results in getting cut more often
	:param data: A numpy array that contains 8 rows of data, each row a numerical equivalent of the amino acid letter
	:param labels: A numpy array that is only 1 dimensional and specifies if the amino acid chain with the same index in the data parameter was cut or not (+1 = cut , -1= not cut)
	
	:return: A scatterplot that plots the portion of chain with a specific amino acid in the x and the percentage of times cut in the Y *NOTE*: No legend included due to there being too many amino acids, any legend would likley be hard to use/identify specific amino acids with
	"""
	
	
	
	unique_aa = list(set(x for l in data for x in l)) ##Making a list of the unique amino acids in the data

	for aa in unique_aa:

		aa_data = [] ##Initializing what will be a list of what the portion of an amino acid each chain is
	
	
		for x in range(len(data)): ##Want to make a function calculates the portion of each chain that is an amino acid
			sub_aa_counter = 0 ##Making a sub amino acid counter

			for y in range(len(data[x])):
				if data[x][y] == aa:
					sub_aa_counter += 1 ##Adding one for every amino acid found in a single chain
				
			aa_data += [sub_aa_counter/8] 

		aa_cut_freq = {} ##Initializng a dicitonary that will store chain portion -> Percentage of how many chains with that portion are cut
	
		for unique in sorted(set(aa_data)): ##Want to initialize dictionaries of each possible portion
			aa_cut_freq[unique] = 0
		
		for i in range(len(aa_data)): ##Recording how many of the portions get cut
	
			if labels[i] == 1: ##Adding to number of chains that have that character only if cleaved
				aa_cut_freq[aa_data[i]] += (1/len(data))*100
			
	
	
		
		plt.scatter(*zip(*aa_cut_freq.items()))
	
	plt.xticks([0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]) ##Tried to use range(0,1,0.125) but won't take floats
	plt.title('Affect of Increasing Portion of Amino Acid in Chain')
	plt.ylabel('Percentage of Chains Cleaved')
	plt.xlabel('Fraction of the Amino Acid Present in Chain')
	plt.show()
	
def learn1(data,labels):
	"""
	Uses SVM machine learning to train a model that can predict if a given amino acid chain would be cut by HIV-1 Protease
	:param data: A numpy array that contains 8 rows of data, each row a numerical equivalent of the amino acid letter
	:param labels: A numpy array that is only 1 dimensional and specifies if the amino acid chain with the same index in the data parameter was cut or not (+1 = cut , -1= not cut)
	
	:return: Accuracy of the model during training and testing. Tests using 50 datapoints from the beginning and 50 from the end of the data (not trained using these)
	"""
	
	##Increasing from 50/-50 to 100/-100 increases testing accuracy
	
	data_training = np.array(data[50:-50].tolist()) ##Training using all data from first 50 to last 50
	data_testing = np.array(data[-50:].tolist() + data[:50].tolist()) ##Testing using first and last 50 obsevations
	labels_training = np.array(labels[50:-50].tolist()) ##Training using all labels/classes from first 50 to last 50
	labels_testing = np.array(labels[-50:].tolist() + labels[:50].tolist()) ##Testing using classes from 50 first and last observations
	
	##Initializing svc and fitting using only training data
	svc = svm.SVC(kernel='linear')
	svc.fit(data_training, labels_training)

	# Compute training accuracy
	correct = 0
	for i in range(len(data_training)):
		prediction = svc.predict([data_training[i]])
		label = labels_training[i]
		if prediction[0] == label:
			correct += 1
	accuracy = correct / len(data_training)

	# Training accuracy
	print('SVC Training accuracy: {:.2f}'.format(accuracy))

	# Compute testing accuracy
	correct = 0
	for i in range(len(data_testing)):
		prediction = svc.predict([data_testing[i]])
		label = labels_testing[i]
		if prediction[0] == label:
			correct += 1
	accuracy = correct / len(data_testing)

	# Testing accuracy
	print('SVC Testing accuracy: {:.2f}'.format(accuracy))		
		
	
def learn2(data,labels):
	"""
	Uses K nearest neighbour machine learning to train a model that can predict if a given amino acid chain would be cut by HIV-1 Protease
	:param data: A numpy array that contains 8 rows of data, each row a numerical equivalent of the amino acid letter
	:param labels: A numpy array that is only 1 dimensional and specifies if the amino acid chain with the same index in the data parameter was cut or not (+1 = cut , -1= not cut)
	
	:return: Accuracy of the model during training and testing. Tests using 50 datapoints from the beginning and 50 from the end of the data (not trained using these)
	"""

		
	data_training = np.array(data[50:-50].tolist()) ##Training using all data from first 50 to last 50
	data_testing = np.array(data[-50:].tolist() + data[:50].tolist()) ##Testing using first and last 50 obsevations
	labels_training = np.array(labels[50:-50].tolist()) ##Training using all labels/classes from first 50 to last 50
	labels_testing = np.array(labels[-50:].tolist() + labels[:50].tolist()) ##Testing using classes from 50 first and last observations
	
	from sklearn.neighbors import KNeighborsClassifier
	
	knn = KNeighborsClassifier()
	knn.fit(data_training, labels_training)

	# Compute training accuracy
	correct = 0
	for i in range(len(data_training)):
		prediction = knn.predict([data_training[i]])
		label = labels_training[i]
		if prediction[0] == label:
			correct += 1
	accuracy = correct / len(data_training)

	# Training accuracy
	print('kNN Training accuracy: {:.2f}'.format(accuracy))

	# Compute testing accuracy
	correct = 0
	for i in range(len(data_testing)):
		prediction = knn.predict([data_testing[i]])
		label = labels_testing[i]
		if prediction[0] == label:
			correct += 1
	accuracy = correct / len(data_testing)

	# Testing accuracy
	print('kNN Testing accuracy: {:.2f}'.format(accuracy))
	
	
## Only used '1625Data.csv' due to it being unsorted. If other datasets that are sorted are used (such as schillingData.csv) testing accuracy of model will be calculated using data that is very similar to data trained with, and accuracy will be higher than actual
	
##TEST CODE, UNCOMMENT TO RUN: 

data,labels = load_data('1625Data.csv') 
viz_1(data, labels)
viz_2(data, labels)
learn1(data,labels)
learn2(data,labels)
