# Code for a study of precursor abundance and transformation rate as a function of distance from the edge and time since deposition

# UPDATES IN PROGRESS
Relevant publication here: 
# System Requirements
MATLAB r2023b and Python 3.11.3 were used in this work. MATLAB is open source and free to use for enrolled students at many educational institutes, or available for purchase. Python is free and open source. Following information details the requirements to run these packages. Python is in general backwards compatible but individual packages may vary in their backwards capibilities. If possible, use versions provided or later, unless major software updates are described in external documentation. 
## Hardware 
All testing in this work was performed on a machine using an M1 processing chip and 16GB available RAM. Compatible, yet untested hardware requirements are listed below.
### MATLAB
RAM: 4GB minimum, 8GB or greater recommended

Processor: 1.8Ghz minimum, faster is preferred. arm support needed.

Disk Space: ~6GB minimum for installation

### Python
RAM: 2GB minimum, 4GB or greater recommended

Processor: any 64-bit arm supported processor compatable with Operating System

Disk Space: ~100MB for installation, not including other required packages

### Total
Roughly 10 GB RAM, using an arm supported processor, and ample disk space > ~6GB is required/recommended to install and use MATLAB and Python (+ packages and dependencies) in this work.
## Operating Systems
All testing in this work was performed using MacOS Sonoma (14). Compatible, yet untested operating systems are listed below.
### MATLAB
MacOS: Big Sur (11) or later

Windows: 10 or later

Linux: Ubuntu 20.04 or later, or other equivalent, supported distributions (CentOS/RHEL 8, Fedora 36, etc.)
### Python
MacOS: Catalina (10) or later

Windows: 10 or later

Linux: Ubuntu 20.04 or later, or other equivalent, supported distributions (CentOS/RHEL 8, Fedora 36, etc.)

# Installation
MATLAB and Python are used independently in this work. MATLAB is used to extract precursor proportions from data processed using Igor. Python is used to calculate the mean distribution of data, perform fits, calculate 1/e lengths, and plot results. The following installation instructions assume MacOS.
## MATLAB

Install instructions are available [here](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html). MacOS machines also require a Java runtime tool, directions also available in preceding link. 

The required tools in MATLAB for this work are

- Image Processing Toolbox
- Optimization Toolbox
- Simulink
- Statistics and Machine Learning Toolbox
- Symbolic Math
- Text Analytics Toolbox

If you have the MATLAB installer with the license (through purchase or institution), you can include the desired toolboxes during the initial installation process by selecting them. Otherwise, after installing MATLAB you can, in the MATLAB interface,  

- Go to the 'Home' tab and click on 'Add-Ons' --> 'Manage Add-Ons' to see all installed toolboxes.
- Go to the 'Home' tab and click on 'Add-Ons' --> 'Get Add-Ons' to search for necessary toolboxes.
Inside of 'Manage Add-Ons' all installed toolboxes are visible.

AAfter all toolboxes are installed, you can run the [demos](https://github.com/zoerechav/Coral_Skeleton_Edge/blob/main/demos/) available for use in the MATLAB interface. 
## Python
Most Linux distributed systems already have Python installed. Check to see if you have Python installed on your local machine using the following command in your terminal:

`python3 --version`

Any version of Python3 or above should be compatible for this work. If you do not have Python on your local machine, install using Homebrew. If you do not have Homebrew, use the following command in terminal (Homebrew is for MacOS users):

`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`

And then install python, which will automatically install the newest version of Python3 (Python2 is outdated and not forward compatable!). If you are not using a MacOS machine and do not have python ( or do not have python and do not like Homebrew ) install directions for python are located [here](https://www.python.org/downloads/).

`brew install python`

Verifying you have python3 and pip3 should mean you are good to go (`python3 --version` and `pip3 --version` in terminal).

Best practice using Python is to create a project and virtual environment where you install all the required packages to perform your function. The packages inside the environment do not change as you update your local machine, thus preserving the software version of the packages you are using.

A project for this work can be called 'coral'. First make the project directory and navigate there.

`mkdir /path/to/coral`

`cd /path/to/coral`

Make your virtual environment

`python3 -m venv venv`

And activate your environment to install the required packages for this work

`source /path/to/coral/venv/bin/activate`

The required packages for this work are located in [versions.py](https://github.com/zoerechav/Coral_Skeleton_Edge/blob/main/versions.py). Depending on the package, exact version specification may not be necessary. To install packages while your virtual environment is activated, use the following commands.

`pip install package (ex: pip install numpy)`
or
`pip install package == version (ex: pip install numpy==1.23.0`

After all packages are installed, you can run the [demos](https://github.com/zoerechav/Coral_Skeleton_Edge/blob/main/demos/) available for use in the terminal. More advanced code editing and visualization is available through the use of an interpreter such as Spyder or Virtual Software (VS) Code.

# Run Demo
## Python 
To run demo scripts in terminal:

`cd /path/to/coral`

`source /path/to/coral/venv/bin/acticate`

Ensure the data folder and scripts are located in /path/to/coral/

`python /path/to/coral/demo_compare_fit.py` and `python /path/to/coral/demo_exponential.py`

