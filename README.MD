

# McGill sea ice model (SIM) code.

## Setup
### clone SIM
connect to chosen machine via aos.meteo.mcgill.ca<br/>
```cd /where/you/want/your/SIM```<br/>
```git clone git@github.com:McGill-sea-ice/SIM.git```<br/>
```cd SIM```<br/>
if needed download forcing:<br/>
```scp <username>@dumbo:/storage2/common/forcing /path/to/forcing``` <br/>
### create link to forcing: 
```ln -s /path/to/forcing forcing``` 

```git branch -a (list branches)```<br/>
```git checkout <choose your branch>```<br/>

### Initial compile:
```./setup_script```<br/>
```export HDF5_DISABLE_VERSION_CHECK=2```<br/>
### Run the code:
```./zoupa < input_file > outputfile &```<br/>
output files for fields are in the output directory<br/>

Subsequent compiles: 
`scons-3`

## PULL REQUESTS
### Testing procedure

Create a `TESTS` directory: 
- `cd /where/you/want/your/TESTS`
- `mkdir TESTS`
- `cd TESTS`

Get the master branch and reference solution in `REF` for reference
- `git clone git@github.com:McGill-sea-ice/SIM.git REF`
- `cd REF`
- `ln -s /path/to/forcing forcing`
- `./testing_script`

Get your own branch in `NEW` and run the model
- `cd ..`
- `git clone git@github.com:McGill-sea-ice/SIM.git NEW`
- `cd NEW`
- `ln -s /path/to/forcing forcing`
- `git checkout <choose your branch>`
- `./testing_script`

Compare the results
- `./test_results`
- copy result of tests (should get 5 PASS).

### steps for PR on github
1. Click on pull request on github.
2. Enter text to describe code modifications and paste result of tests.
3. Assign someone for reviewing the code.
4. Assign lbtremblay (and JFLemieux73 for a little bit) for approving the pull request.
