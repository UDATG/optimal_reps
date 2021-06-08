# optimal_reps
Optimal cycle representatives

## Installing CBC
`brew tap coin-or-tools/coinorbrew`

`install coin-or-tools/coinor/cbc`

https://github.com/coin-or/Cbc

## Installing Gurobi


### Mac

* make sure you do not move the gurobi license file from the place it's first saved to
* you'll need to set some "environmental variables"
	
	1. we will do this by modifying a `.bash_profile` folder 
		* there are several other files you can modify to update an environmental variable, each with different effects (a few examples of such files include `/.bashrc`, `/etc/profile`, `~/.profile`, `~/.zprofile`)
			* here are some partially overlapping online resources:
				* for Mac OS  [here](https://youngstone89.medium.com/setting-up-environment-variables-in-mac-os-28e5941c771c) 
				* for Mac, Linux, and Windows: [a nice discussion](https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7), 
				* for `bash`, `zsh`, or `tcsh`  [a nice stack exchange](ttps://unix.stackexchange.com/questions/21598/how-do-i-set-a-user-environment-variable-permanently-not-session)

    2. check to see if you have a file `/.bash_profile`
        a. You may have hit "Command + Shift + ." to show hidden files
        b. If you don't already have a `.bash_profile` file, then create one with a code editor like atom, sublime, vim, etc.
    3. add the following lines to your `.bash_profile` file:
    
        ```
        # Environmental variables for Gurobi
        export GRB_LICENSE_FILE="path/to/gurobi.lic" # the license file
        export GUROBI_HOME="path/to/installation/folder" 
        ```
    4. open a terminal and run
    
        ```
        source ~/.bash_profile
        ```
       to refresh your environmental variables with the new addition
    5. run
    
       ```
       printenv
       ```
       to double check that the new variables have been added

#### Set System Environment Variable
1. Build a plist file (https://discussions.apple.com/thread/7814747)
2. Setup, activate and launch the plist as instructed by
      https://support.shotgunsoftware.com/hc/en-us/articles/219042108-Setting-global-environment-variables-on-OS-X
      
* The set_objective function in Gurobi doesn't seem to work. However, the function can compile without this line.

## Examples of homology computations

* There are example files in `src/bin`
* Suppose you want to compile and run the file `demo_rational.rs`, which lives in `src/bin`: 

    ```
    cd path/to/crate
    cargo build --bin demo_rational 
    ./target/debug/demo_rational
    ```
* Alternatively, you can do

    ```
    cd path/to/crate
    cargo run --bin demo_rational 
    ```
* If you want to compile/run a new example, say `demo_z5.rs`, then save `demo_z5.rs` to `src/bin` and follow the steps above.
