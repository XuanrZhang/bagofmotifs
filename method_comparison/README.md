In this file folder, we provided all code for method comparison, including three neural network methods.

- Basset-2016 (see '')
- DeepSTARR-2020 (see '')
- DeepMEL-2020 (see '')

For a fair comparion, we used 500bp sequences from the center of peaks and their reverse_complement sequences to double our samples. 
For each method, we reimplemented their methods according to their works, and only modified the output layer (n=17) to predict 17 cell-type specific enhancers.
More details of those methods can be found in their github.

In this work, we trained model on mouse development data (E8.25) with arguments (Epochs = 100; Batchsize = 128; earlystopping =10), more details of model training and evaluation can be found in the scirpts.

If you don't want to re-train the models, you can directly re-load our trained model files (see ./models file folder)


Any problems during your running, please email to x.zhang@victorchang.edu.au

