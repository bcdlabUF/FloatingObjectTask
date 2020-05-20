# FloatingObjectTask
scripts to process eeg data from the floating objects task

Stimuli_FOT.m is the script used to present the stimuli during testing.

If you do ICA, FOT_Step1_filter_ICA.m filters the data and performs ICA. Then use FOT_Step2_split_clean.m epochs the data and saves a separate file for each condition. Then each new file is cleaned and an average reference is computed. Step 2 uses FOT_Split_Cond_batch.m and BEES_clean_data.m

If you don't use ICA, then use Step1_FOT_pipeline.m which filters the data, epochs, saves a separate file, cleans the data, and references.

Step2_Cohere_and_Amp.m segments each trial into 1s segments, removes noisy segments, computes coherence and amplitude. This pipeline uses cohere_baby.m.
