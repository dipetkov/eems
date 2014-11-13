
# EEMS: python-pipeline

this is a basic pipeline to run eems. It is writen as an ipython notebook, with
the potential goal of making the resulting python script directly runnable.

Right now, the goal is to make a working prototype that does the following
steps:

 1. reads VCF
    - also read other files
 2. calculate distance matrix
 3. run eems
 4. create figures



###Install:

I am still experimenting with modules, but to install the notebook, simply

```
sudo apt-get install ipython-notebook python-numpy python-matplotlib python-pandas
ipython notebook
``` 
then  view the notebook in a webbrowser
