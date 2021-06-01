# optimal_reps
Optimal cycle representatives

## Installing CBC
`brew tap coin-or-tools/coinorbrew`

`install coin-or-tools/coinor/cbc`

https://github.com/coin-or/Cbc

## Installing Gurobi


### Mac

* make sure you do not move the gurobi license file from the place it's first saved to
* you'll need to set some "environmental variables", and you'll do this by modifying a `.bash_profile` folder
    1. check to see if you have a file `/.bash_profile`
        a. You may have hit "Command + Shift + ." to show hidden files
        b. If you don't already have a `.bash_profile` file, then create one with a code editor like atom, sublime, vim, etc.
    2. add the following lines to your `.bash_profile` file:
    
        ```
        # Environmental variables for Gurobi
        export GRB_LICENSE_FILE="path/to/gurobi.lic" # the license file
        export GUROBI_HOME="path/to/installation/folder" 
        ```
    3. open a terminal and run
    
        ```
        source ~/.bash_profile
        ```
       to refresh your environmental variables with the new addition
    4. run
    
       ```
       printenv
       ```
       to double check that the new variables have been added


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
