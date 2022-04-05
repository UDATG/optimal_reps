use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind,FactoredComplexBlockCsm};
use exhact::clique::Simplex;
use num::rational::Ratio;
use std;
extern crate gurobi;
use gurobi::*;
use sprs::{CsMat, CsVec};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::hash::Hash;
use optimal_representatives::simplex_bar::{simplex_barcode};
use ndarray::Array;
use ndarray::array;
use ndarray_npy::write_npy;
use num::integer::Integer;
type Coefficient = Ratio<i16>;
use ordered_float::OrderedFloat;

fn ordered_floats( v : Vec<f64> ) -> Vec< OrderedFloat<f64> > {
    let u : Vec<_> = v.into_iter().map(OrderedFloat).collect(); 
    return u
}

fn ordered_floats_nested(v: Vec<Vec<f64>>) -> Vec< Vec< OrderedFloat<f64> > > {
    return v.into_iter().map( ordered_floats ).collect();
}
fn main() {  
    let mut f = BufReader::new(File::open("data_text\\3_d_dis_mat.txt").unwrap());
    let mut s = String::new();

     // for the input as Vec-of-Vec square symmetric matrix
     let arr: Vec<Vec<f64>> = f.lines()
     .map(|l| l.unwrap().split(char::is_whitespace)
          .map(|number| number.parse().unwrap())
          .collect())
     .collect();
   

    let dismat = ordered_floats_nested(arr.clone());
    
    // set the max dimension to compute persistent homology
    let mut dim_input = String::new();
    println!("Set the dimension of the persistent homology:");
    let mut dim = 1;
    std::io::stdin().read_line(&mut dim_input);
    dim = dim_input.trim().parse::<usize>().unwrap();


    // set the maximum dissimilarity threshold
    let mut max_threshold = f64::INFINITY;
    for row in arr.clone(){
        let row_max = row.iter().cloned().fold(0./0., f64::max);
        if row_max < max_threshold{
            max_threshold = row_max;
        }
    }
    
    let maxdis: OrderedFloat<f64> = OrderedFloat(max_threshold);

    // create a "ring object" representing the field of rational numbers
    let ringmetadata = exhact::matrix::RingMetadata{
    ringspec: RingSpec::Rational,
    identity_additive: Ratio::new(0, 1),
    identity_multiplicative: Ratio::new(1, 1)
    };

    
    // build and factor the filtered chain complex
    let chx = exhact::clique::CliqueComplex {
        // the distance/dissimilarity matrix
        dissimilarity_matrix: dismat.clone(), 
        // threshold to stop the filtration
        dissimilarity_value_max: maxdis, 
        // sets "safeguards" on dimension; we'll get warnings if we try to 
        // get boundary matrices in dimension higher than dim+1
        safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(), 
        // set the default major dimension (for sparse matrices) to be row
        major_dimension: MajorDimension::Row, 
        // indicates we want Z/3Z coefficients
        ringmetadata: ringmetadata, 
        // don't worry about this
        simplex_count: Vec::new() 
    };
    
    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);

    let mut birth_vec: Vec<f64> = Vec::new();
    let mut death_vec: Vec<f64> = Vec::new();

    // obtain a list of (birth_edge, death_triangle) pairs for the nonzero bars 
    let simplex_bar = simplex_barcode( &factored_complex, dim );
    

    for j in 0..simplex_bar.len(){
        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        
        birth_vec.push(chx.key_2_filtration(&birth).into());
        death_vec.push(chx.key_2_filtration(&death).into());
    }

    let birth_arr = Array::from_vec(birth_vec);
    let death_arr = Array::from_vec(death_vec);
    

    write_npy("simplex_bar/3_d_1_birth_time.npy", &birth_arr);
    write_npy("simplex_bar/3_d_1_death_time.npy", &death_arr);
}