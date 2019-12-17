#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################
##### ALFRESCO ff-DNN builder #####
###################################

"""
This script grabs pre-constructed data from a standard ALFRESCO directory environment; it can be run when the datasets outlined in the 
dictionary 'dataset_paths' are present. The script will, per dataset, construct simple feed-forward DNNs learned on 
ddGoffset values present in the input files using SKLearn to preprocess and tensorflow as a learning and distribution platform.

A bayesian optimiser (skopt) loops (n=40-50, ideally) through each replicate to search hyperparameter space, optimising on layers, neurons, 
batch size and the Adam parameters beta1, beta2 and epsilon. In each optimisation process the best-performing model will be saved 
to 'opt_output/{ensemble_ID}_{split}_ALFRESCO_TopPerform_architecture.json' and 'opt_output/{ensemble_ID}_{split}_ALFRESCO_TopPerform_weights.h5',
from where it can be read in using other scripts to predict on external sets the ALFRESCO predictor script leverages the saved models as an ensemble. 
Note that using skopt any set of hyperparameters can be used to optimise on. For each optimal model, internal and external validation results 
are written to 'opt_output/{split}_TopPerformer_internalVal_df.csv' and 
'opt_output/{split}_TopPerformer_externalVal_df.csv', respectively; these can be plotted using external scripts.

Finally, the script also outputs the cumulative minima of prediction error (on the validation set) per call and replica for each dataset.
Running plot_convergence.py draws a pretty graph (seaborn) of (potential) convergence in hyperparameter space.

Due to skopt's acquisition function scaling issues, using more than 50-60 calls can become exponentially expensive and slow down the program significantly.
We recommend using 40 calls and check if convergence has been reached by using the plot_convergence.py script.

"""

# TF-related imports & some settings to reduce TF verbosity:
import os
os.environ["CUDA_VISIBLE_DEVICES"]="1, 2, 3"	# current workstation contains 4 GPUs; exclude 1st
import tensorflow as tf 
from tensorflow import keras
tf.logging.set_verbosity(tf.logging.ERROR)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# SciKit-Optimize:
import skopt
from skopt import gp_minimize, forest_minimize
from skopt.space import Real, Categorical, Integer
from skopt.plots import plot_convergence
from skopt.plots import plot_objective, plot_evaluations
from tensorflow.python.keras import backend as K
from skopt.utils import use_named_args

# General imports:
import glob
import shutil
import subprocess
import numpy as np
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import seaborn as sns
import matplotlib.cbook
import time
start_time = "Start: "+str(time.ctime())
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# Misc. imports:
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
from sklearn.metrics import mean_absolute_error
from sklearn.decomposition import PCA
from scipy import stats
import statistics
import pickle
print(
"\n"
"\n"
"##############################################################\n"
"######################  ALFRESCO - Build #####################\n"
"##############################################################\n"
"###### Target-based-Cross-Validation feed-forward-DNN ########\n"
"# model generator using Bayesian hyperparameter optimisation #\n"
"##############################################################\n")
print("Imports complete.\n")

# set dataset to optimise on:
dataset_paths = {
 '../datasets/trainingsets_compiled/dataset_1.csv': "1",
'../datasets/trainingsets_compiled/dataset_2.csv': "2",
'../datasets/trainingsets_compiled/dataset_3.csv': "3",
'../datasets/trainingsets_compiled/dataset_12.csv': "12",
'../datasets/trainingsets_compiled/dataset_13.csv': "13",
'../datasets/trainingsets_compiled/dataset_23.csv': "23",
'../datasets/trainingsets_compiled/dataset_123.csv': "123",
'../datasets/trainingsets_compiled/dataset_noise.csv': "Noise"
}

# list of ligand names in compiled dataset files:
# update when new targets are added!!

targets = [
			"HSP90",
			"jnk1",
			"FXR", 
			"ACK1",
			"throm_jm", 
			"tyk", 
			"throm_schr",
			#"test", 
			"BACE", 		
			"cdk2",
			"mcl1",
			"ptp1b"
			]

# set data processing configurations:
PCA_threshold = 0.95				# Keeps n dimensions for x variance explained
replicates =30						# Number of replicates per subject model
n_calls = 40						# Number of Bayesian optimisation loops for hyperparameter optimisation, 40 is best for convergence, > 60 scales to very expensive
startpoint_BO = np.inf				# Point to consider top-performing model from (MAE/MAD); 1.0 = no improvement on test-set variance
ensemble_size = 10					# Amount of top-scoring models to retain per fold-dataset combination
runtime_estimation = replicates * 37.6	* len(dataset_paths)/60			

print(
	"Program is set with the following settings:\n"
	"Feature sets:", str(len(dataset_paths)),"\n"
	"Protein (i.e. perturbation) sets:", str(len(targets)),"\n"
	"PCA variance retained threshold = "+str(PCA_threshold)+"\n"			
	"Replicates = "+str(replicates)+"\n"			
	"Model ensemble size = "+str(ensemble_size)+"\n"
	"\n"
	)
if n_calls == 40:
	print("Estimated runtime:", str(round(runtime_estimation, 1)), "hours")
else:
	print("Could not estimate runtime because n_calls is not 40.")

def TranslateTargetNames(input_name):
	newname = input_name
	
	if "throm_jm" in newname:
		newname = "THROMBIN-JM"
	if "throm_schr" in newname:
		newname = "THROMBIN-SCHR"
	if "cdk2" in newname:
		newname = "CDK2"
	if "FXR" in newname:
		newname = "FXR"			
	if "ACK1" in newname:
		newname = "ACK1"
	if "tyk" in newname:
		newname = "TYK2"
	if "jnk1" in newname:
		newname = "JNK1"
	if "HSP90" in newname:
		newname = "HSP90"
	if "BACE" in newname:
		newname = "BACE"
	if "mcl1" in newname:
		newname = "MCL1"
	if "ptp1b" in newname:
		newname = "PTP1B"


	if "test" in newname:
		newname = "2ndBACE"
		
			
	return newname

def NormaliseDatasets(collection):
	# process input nested lists of datasets (per target):
	train_dataset = collection
	
	print("Normalising..")
	# Calculate statistics, compute Z-scores, clean:
	train_stats = train_dataset.describe()

	train_stats.pop("ddG_offset")
	train_stats = train_stats.transpose()

	train_labels = train_dataset.pop('ddG_offset')


	def norm(x):
		return (x - train_stats['mean']) / train_stats['std']

	# Normalise and return seperately:
	normed_train_data = norm(train_dataset).fillna(0).replace([np.inf, -np.inf], 0.0)
	
	
	return [normed_train_data, train_labels]



def ReduceFeatures(normalised_collection, PCA_threshold, split):
	print("Computing PCA, reducing features up to "+str(round(PCA_threshold*100, 5))+"% VE..")
	training_data = normalised_collection
	
	# Initialise PCA object, keep components up to x% variance explained:
	PCA.__init__
	pca = PCA(n_components=PCA_threshold)


	# Fit to and transform training set:			
	train_postPCA = pd.DataFrame(pca.fit_transform(training_data))

	print("# of PCA features after reduction: "+str(len(train_postPCA.columns)))

	train_postPCA.index = training_data.index
	# pickle pca object to file so that external test sets can be transformed accordingly (see https://stackoverflow.com/questions/42494084/saving-large-data-set-pca-on-disk-for-later-use-with-limited-disc-space)
	pickle.dump(pca, open("./opt_output/pca_featureset_"+str(split)+".p", "wb"))
	return train_postPCA	# return list with test_postPCA when needed


def SplitDatasets(dataset, targets, labels):
	print("Splitting data per target..")


	dataset = pd.concat([dataset, labels], axis=1)

	# split datasets by regex, store into dict:
	target_splits_dict = {}
	
	for target_name in targets:
		split_set = dataset[dataset.index.str.contains(target_name, regex=False)]
		target_splits_dict[target_name] = split_set

	# construct train-test combinations for all target datasets:
	CV_splits = []
	for key, test_set in target_splits_dict.items():
		training_set = [ rows for target, rows in target_splits_dict.items() if key != target ]
		training_set = pd.concat(training_set)


		train_labels = training_set["ddG_offset"]
		test_labels = test_set["ddG_offset"]
		
		training_set = training_set.drop("ddG_offset", axis=1)
		test_set = test_set.drop("ddG_offset", axis=1)

		CV_splits.append([[training_set, test_set], [train_labels, test_labels]])

	print("Done. Initialising tb-CV optimisation loops..")
	return CV_splits


def FF_DNN_KERAS(dataframe, dataset_name, iteration):
	tf.keras.backend.clear_session()
	tf.reset_default_graph()
	model_bucket = []
	# Display training progress by printing a single dot per epoch:
	class PrintDot(keras.callbacks.Callback):
	  def on_epoch_end(self, epoch, logs):
	    if epoch % 100 == 0: print('')
	    print('.', end='')

	# Set early stopping variable:
	early_stopping = keras.callbacks.EarlyStopping(
													monitor='val_loss', 
													mode='min', 
													patience=20,
													verbose=0)



	# Retrieve datasets:
	train_postPCA, test_postPCA, train_labels, test_labels = dataframe[0][0], dataframe[0][1], dataframe[1][0], dataframe[1][1]

	

	# Build keras DNN using global params:
	def create_model(
		num_dense_layers_base, 
		num_dense_nodes_base,
		num_dense_layers_end, 
		num_dense_nodes_end, 
		activation,
		adam_b1,
		adam_b2,
		adam_eps,
		num_batch_size):


		model = keras.Sequential()

	# Add input layer of length of the dataset columns:
		model.add(keras.layers.Dense(len(train_postPCA.columns), input_shape=[len(train_postPCA.keys())]))

	# Generate n number of hidden layers (base, i.e. first layers):
		for i in range(num_dense_layers_base):
			model.add(keras.layers.Dense(num_dense_nodes_base,
			activation=activation
			))

	# Generate n number of hidden layers (end, i.e. last layers):
		for i in range(num_dense_layers_end):
			model.add(keras.layers.Dense(num_dense_nodes_end,
			activation=activation
			))

	# Add output layer:

		model.add(keras.layers.Dense(1, activation=keras.activations.linear))

		optimizer = tf.keras.optimizers.Adam(lr=0.0001, beta_1=adam_b1, beta_2=adam_b2, epsilon=adam_eps)

		model.compile(
			loss="mae",
			optimizer=optimizer,
			metrics=["mae"]
			)
		return model


	# Set hyperparameter ranges, append to list:
	dim_num_dense_layers_base = Integer(low=1, high=2, name='num_dense_layers_base')
	dim_num_dense_nodes_base = Categorical(categories=list(np.linspace(5,261, 10, dtype=int)), name='num_dense_nodes_base')
	dim_num_dense_layers_end = Integer(low=1, high=2, name='num_dense_layers_end')
	dim_num_dense_nodes_end = Categorical(categories=list(np.linspace(5,261, 10, dtype=int)), name='num_dense_nodes_end')


	#dim_activation = Categorical(categories=[tf.keras.activations.relu], name='activation')
	dim_adam_b1 = Categorical(categories=list(np.linspace(0.8,0.99,11)), name="adam_b1")
	dim_adam_b2 = Categorical(categories=list(np.linspace(0.8,0.99,11)), name="adam_b2")
	dim_adam_eps = Categorical(categories=list(np.linspace(0.0001, 0.5, 11)), name="adam_eps")
	dim_num_batch_size = Categorical(categories=list(np.linspace(32, 128, 7, dtype=int)), name='num_batch_size')
	
	dimensions = [
				dim_num_dense_layers_base,
				dim_num_dense_nodes_base,
				dim_num_dense_layers_end,
				dim_num_dense_nodes_end,
				dim_adam_b1,
				dim_adam_b2,
				dim_adam_eps,
				dim_num_batch_size]	
	
	@use_named_args(dimensions=dimensions)
	def fitness(
		num_dense_layers_base, 
		num_dense_nodes_base, 
 
		num_dense_layers_end, 
		num_dense_nodes_end,
		adam_b1,
		adam_b2,
		adam_eps,
		num_batch_size):
    # Create the neural network with these hyper-parameters:
		model = create_model(
							num_dense_layers_base=num_dense_layers_base,
							num_dense_nodes_base=num_dense_nodes_base,
							num_dense_layers_end=num_dense_layers_end,
							num_dense_nodes_end=num_dense_nodes_end,
							activation=tf.keras.activations.relu,
							adam_b1=adam_b1,
							adam_b2=adam_b2,
							adam_eps=adam_eps,
							num_batch_size=num_batch_size)


		#print("Fitting model..")
		history = model.fit(
			train_postPCA, train_labels,
		epochs= 1000, 
		validation_data = (test_postPCA, test_labels),
		verbose=0,
		callbacks=[
					early_stopping, 
					#PrintDot(),			# uncomment for verbosity on epochs
					], 		
		batch_size=121)

		hist = pd.DataFrame(history.history)
		hist['epoch'] = history.epoch
		MAE = hist["val_mean_absolute_error"].tail(10).mean()
		MAD_testset = test_labels.mad()

		MAEMAD = MAE/MAD_testset
		print("MAE/MAD:",MAEMAD)

		# calculate some statistics on test set:
		prediction = model.predict(test_postPCA)
		prediction_list = [ item[0] for item in prediction ]

		perts_list = test_labels.index.tolist()
		exp_list = test_labels.values.tolist()

		slope, intercept, r_value, p_value, std_err = stats.linregress(prediction_list, exp_list)
		tau, p_value = stats.kendalltau(prediction_list, exp_list)


		# For plotting test set correlations:
		tuples_result = list(zip(perts_list, exp_list, prediction_list))
		nested_list_result = [ list(elem) for elem in tuples_result ]

	# Append data with best performing model.
	# Data contains the MAE/MAD score, protein target, iteration,
	# tau, r value, the keras DNN model, the internal validation plot 
	# and the data for external validation:
		global startpoint_MAEMAD

		if MAEMAD < startpoint_MAEMAD:
			startpoint_MAEMAD = MAEMAD
			model_bucket.append([MAEMAD, dataset_name, iteration, tau, r_value, model, hist, nested_list_result])

			# # write all model files:
			# # Slightly hacky but TF's backend voids model parameters when the model is saved as a variable
			# # in order to retain the top performing model. From these temporary model files, all but the 
			# # top-performing model will be deleted from the system at the end of this script.
			if not os.path.exists("./opt_tmp"):
				os.makedirs("./opt_tmp")

			model.save_weights("opt_tmp/"+str(iteration)+"_"+str(translated_subject)+"_ALFRESCO_TopPerform_weights.h5")
			with open("opt_tmp/"+str(iteration)+"_"+str(translated_subject)+"_ALFRESCO_TopPerform_architecture.json", "w") as file:
				file.write(model.to_json())
		
		del model
		tf.keras.backend.clear_session()
		K.clear_session()		
		
		return MAEMAD

	# Bayesian Optimisation to search through hyperparameter space. 
	# Prior parameters were found by manual search and preliminary optimisation loops. 
	# For running just dataset 13x500 calls, optimal hyperparameters from 150 calls were used as prior.
	default_parameters = [2, 33, 1, 90, 0.971, 0.895, 1.0000e-04, 112]
	print("###########################################")
	print("Created model, optimising hyperparameters..")
	
	search_result= gp_minimize(func=fitness,
								dimensions=dimensions,
								acq_func='EI', #Expected Improvement.
								n_calls=n_calls,
								x0=default_parameters)


	print("###########################################")
	print("Concluded optimal hyperparameters:")
	print(search_result.x)

	print("###########################################")

	# return skopt object and highest scoring model for this replicate:
	return search_result, model_bucket[-1]



##########################################################
##########################################################
######												######
######				Function calls					######
######												######
##########################################################
##########################################################


# initiate empty DF to fill with cumulative minima 
cumulative_MAEs = pd.DataFrame()
cumulative_MAEtauR_CV = pd.DataFrame()
# clean slate opt_output:
if os.path.exists("./opt_output"):
	subprocess.call('rm ./opt_output/*', shell=True)
if not os.path.exists("./opt_output"):
	os.mkdir("./opt_output")


# initiate log file:
with open("opt_output/logfile.txt", "w") as file:
					writer = csv.writer(file, delimiter='\t')
					writer.writerow(["###########Starting tb-CV BO.###########"])
					writer.writerow(["PCA threshold: "+str(PCA_threshold)])
					writer.writerow(["n replicates: "+str(replicates)])
					writer.writerow(["n models in ensemble: "+str(ensemble_size)])
					writer.writerow(["n calls (BO): "+str(n_calls)])
					writer.writerow(["Started program at: "+time.ctime()])

# loop over input file dict:
for dataset, split in dataset_paths.items():
	bucket_df = pd.DataFrame()
	mae_results_per_fold = [["Subject", "MAE", "Replicate"]]
	MAEtauR_results_per_fold = [["Subject", "Correlation Coefficient", "Dataset", "Correlation metric"]]

	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	print(time.ctime())
	print("Working on dataset: "+split)
# construct raw dataset, make sure all values are floats; shuffle rows:
	collection = pd.read_csv(dataset, index_col="Perturbation")
	collection = collection.apply(pd.to_numeric).astype(float).sample(frac=1)

	normalised_dataset, labels = NormaliseDatasets(collection)
	
# reduce features (and sparsity) using PCA:
	preprocessed_data = ReduceFeatures(normalised_dataset, PCA_threshold, split)
	
# Split dataset per target:
	split_datasets = SplitDatasets(preprocessed_data, targets, labels)

# start loop on each protein fold:
	for cross_val_split in split_datasets:
		
		
		subject = cross_val_split[0][1].index[0].split(">")[0]

		translated_subject = TranslateTargetNames(subject)
		print("Working on subject:", translated_subject)
		

	# Collect MAEs for statistics:
		MAEs_per_split = []
		models_per_replicate = []
		for i in range(replicates):
		# run tb-CV:
			# reset MAEMAD startpoint per replicate:
			startpoint_MAEMAD = startpoint_BO
			OptimizeResult, top_performer = FF_DNN_KERAS(cross_val_split, split, i)

			models_per_replicate.append(top_performer)

		# construct, cummin and concatenate results of this replicate to the other replicates in the loop:
			split_columns = { 
							"Dataset" : str(split), 
							"MAE/MAD" : OptimizeResult.func_vals,
							"Subject": translated_subject}
			result_df = pd.DataFrame(split_columns).cummin()
			bucket_df = pd.concat([bucket_df, result_df])
		# tag data with the dataset type (i.e. descriptor set), add to complete results:
			bucket_df["Dataset"] = str(split)
			cumulative_MAEs = pd.concat([cumulative_MAEs, bucket_df])


		# retrieve statistics for this replicate:					
			tau = top_performer[3]
			r_value = top_performer[4]
			MAE = top_performer[0]

			MAEtauR_results_per_fold.append([translated_subject, r_value, split, "Pearson's-r"])
			MAEtauR_results_per_fold.append([translated_subject, tau, split, "Kendall's-tau"])
			MAEtauR_results_per_fold.append([translated_subject, MAE, split, "MAE/MAD"])
		
		# write update to log file:
			with open("opt_output/logfile.txt", "a") as file:
				writer = csv.writer(file, delimiter='\t')
				writer.writerow(["Finished "+translated_subject+", dataset "+split+", replicate "+str(i+1)+" at "+str(time.ctime())])

		# make ensemble of best models; pick n replicates' top performing models:
		
		models_per_replicate = sorted(models_per_replicate, key=lambda x: x[0])

		ensemble_collection = models_per_replicate[:ensemble_size]
		
		i=1
		for best_model_collection in ensemble_collection:
			
			opt_replicate = str(best_model_collection[2])
			nested_list_result_internal = best_model_collection[6]

			nested_list_result_external = best_model_collection[7]


	# For this best model, retrieve model files, plot internal validation and write external validation to file:
			if not os.path.exists("./opt_output"):
				os.mkdir("./opt_output")	
				
			# with the known optimal replicate #, isolate model files from opt_tmp and move to opt_output:
			# rename so that name contains name of the feature set instead of the replicate:
			os.rename(
				"opt_tmp/"+opt_replicate+"_"+translated_subject+"_ALFRESCO_TopPerform_architecture.json",
				"opt_output/model"+str(i)+"_"+split+"_"+translated_subject+"_ALFRESCO_TopPerform_architecture.json"
				)
			os.rename(
				"opt_tmp/"+opt_replicate+"_"+translated_subject+"_ALFRESCO_TopPerform_weights.h5",
				"opt_output/model"+str(i)+"_"+split+"_"+translated_subject+"_ALFRESCO_TopPerform_weights.h5"
				)
			i+=1
		# to keep things clean, remove ./opt_tmp:
		shutil.rmtree("./opt_tmp/")

		# write internal validation DF: 
		nested_list_result_internal.to_csv("opt_output/"+str(split)+"_"+str(translated_subject)+"_TopPerformer_internalVal_df.csv")

		# write external validation DF:
		with open("opt_output/"+str(split)+"_"+str(translated_subject)+"_TopPerformer_externalVal_df.csv", "w") as file:
			writer = csv.writer(file)
			writer.writerow(["Perturbation", "Experimental ddGoffset (kcal/mol)", "Predicted ddGoffset (kcal/mol)", "Subject"])
			for row in nested_list_result_external:
				writer.writerow(row + [translated_subject])


	
	MAEs_CV = pd.DataFrame(mae_results_per_fold[1:], columns=mae_results_per_fold[0])
	MAEtauR_CV = pd.DataFrame(MAEtauR_results_per_fold[1:], columns=MAEtauR_results_per_fold[0])
	cumulative_MAEtauR_CV = pd.concat([cumulative_MAEtauR_CV, MAEtauR_CV])


cumulative_MAEtauR_CV.to_csv("output/tb-CV_MAEtauR_outputs_50x40.csv", index=False)
cumulative_MAEs.to_csv("output/tbCV_BO_MAE_50x40.csv")
print("Success, wrote all files to opt_output/.")



