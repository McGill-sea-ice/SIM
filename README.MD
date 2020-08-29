

McGill sea ice model (SIM) code.

connect to dumbo or other machine via zephyr<br/>
cd where you want to go
git clone git@github.com:McGill-sea-ice/SIM.git
cd SIM
git branch -a (list branches)
git checkout choose your branch
./setup_script (to compile and create link for forcing) - ignore warnings
./zoupa < input_file > outputfile &
output files for fields are in the output directory


******* PULL REQUESTS *******

-Run 1 day exp with new code and 1 day exp with reference code.
-Run check_BFB.py routine to verify if it is bit for bit (BFB).
-Click on pull request.
-Enter text to describe code modifications and mention if it is BFB or not. 
-Assign someone for reviewing the code.
-Assign JFLemieux73 for approving the pull request. 
