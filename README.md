# optimal_reps
Optimal cycle representatives

## Installing CBC
`brew tap coin-or-tools/coinorbrew`
`install coin-or-tools/coinor/cbc`
https://github.com/coin-or/Cbc

## Examples

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
