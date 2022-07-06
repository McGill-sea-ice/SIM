

McGill sea ice model (SIM) code.

connect to dumbo or other machine via zephyr<br/>
cd where you want to go<br/>
git clone git@github.com:McGill-sea-ice/SIM.git<br/>
cd SIM<br/>
git branch -a (list branches)<br/>
git checkout \<choose your branch\><br/>
. setup_script (to compile and create link for forcing)<br/>
export HDF5_DISABLE_VERSION_CHECK=2<br/>
./zoupa < input_file > outputfile &<br/>
output files for fields are in the output directory<br/>


******* PULL REQUESTS *******<br/>

mkdir TESTS<br/>
cd TESTS<br/>
git clone git@github.com:McGill-sea-ice/SIM.git REF <br/>
cd REF<br/>
. testing_script<br/>
cd ..<br/>
git clone git@github.com:McGill-sea-ice/SIM.git NEW <br/>
cd NEW<br/>
git checkout \<choose your branch\><br/>
. testing_script<br/>
./test_results<br/>
copy result of tests (should get 5 PASS).<br/>
Click on pull request on github.<br/>
Enter text to describe code modifications and paste result of tests.<br/>
Assign someone for reviewing the code.<br/>
Assign lbtremblay (and JFLemieux73 for a little bit) for approving the pull request.<br/>
