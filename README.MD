

McGill sea ice model (SIM) code.

connect to dumbo or other machine via zephyr<br/>
cd where you want to go<br/>
git clone git@github.com:McGill-sea-ice/SIM.git<br/>
cd SIM<br/>
git branch -a (list branches)<br/>
git checkout choose your branch<br/>
./setup_script (to compile and create link for forcing)<br/>
export HDF5_DISABLE_VERSION_CHECK=2<br/>
./zoupa < input_file > outputfile &<br/>
output files for fields are in the output directory<br/>


******* PULL REQUESTS *******<br/>

-Run 1 day exp with new code and 1 day exp with reference code.<br/>
-Run check_BFB.py routine to verify if it is bit for bit (BFB).<br/>
-Click on pull request.<br/>
-Enter text to describe code modifications and mention if it is BFB or not.<br/>
-Assign someone for reviewing the code.<br/>
-Assign JFLemieux73 for approving the pull request.<br/>
