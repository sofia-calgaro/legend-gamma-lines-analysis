This code works fine with version v1.0.0 of BAT. Follow the following steps to install BAT and run the script correctly. The use of the container is left out since there are computational resource restrictions in it that do not allow to perform high precision Bayesian fits.

# How to install BAT v1.0.0
Go to the [BAT github repository ](https://github.com/bat/bat) and download the main branch (as a ZIP \[then in case use `tar -xvf bat_code.tar.gz`\] or a git clone). 
Once you have your BAT folder, do the following (put your username where necessary):

``` bash
$ cd bat_downloaded_folder
$ ./autogen.sh 
$ ./configure --prefix=/home/your_username/bat_install
$ make
$ make install
```

BAT is installed if you see on terminal the message
``` bash
============================================================
                 BAT INSTALLATION SUCCESSFUL!     
```

You'll see the message also suggest you to set up your environment:
``` bash
In the bash, you can set up your environment with

export PATH="/home/calgaro/bat_install/bin:$PATH"
export LD_LIBRARY_PATH="/home/calgaro/bat_install/lib:$LD_LIBRARY_PATH"
export CPATH="/home/calgaro/bat_install/include:$CPATH"
export PKG_CONFIG_PATH="/home/calgaro/bat_install/lib/pkgconfig:$PKG_CONFIG_PATH"
```
Which means, you jsut have to open `~/.bashrc` and copy paste the 4 above lines.

# How to run a gamma fit
Once BAT is installed, you are ready to run the code!

To do that, you need a Makefile. Make sure the paths to source and libraries are correct. 

In particular, change the following lines with your username:
``` bash
BAT_INCLUDE = /home/your_username/bat_install/include 
BAT_LIB = /home/your_username/bat_install/lib
```

ROOT files should be located in the same place for all users. In case, modify the following lines:
``` bash
CXXFLAGS += -I/lfs/l1/legend/software/root/v06.14/include
```

Once done, run on terminal
``` bash
make
```
to generate an executable, you can later run as
``` bash
./runGammaAnalysis
```

Everything contained in the main function of "runGammaAnalysis.cxx" will be executed!
