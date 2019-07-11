# FloatingObjectTask
scripts to process eeg data from the floating objects task

Step1_FOT_pipeline.m filters the data, epochs the data, and saves a separate file for each condition. Then each new file is cleaned and an average reference is computed. The pipeline uses FOT_Split_Cond_batch.m and BEES_clean_data.m

Step2_Cohere_and_Amp.m segments each trial into 1s segments, removes noisy segments, computes coherence and amplitude. This pipeline uses cohere_baby.m.
