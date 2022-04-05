use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind,FactoredComplexBlockCsm};
use exhact::clique::Simplex;
use num::rational::Ratio;
use std;
use std::path::Path;
use std::fs;
extern crate gurobi;
use gurobi::*;
use sprs::{CsMat, CsVec};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::hash::Hash;
use optimal_representatives::simplex_bar::{simplex_barcode};
use exhact::solver::multiply_hash_smoracle; // multiply a hashmap by a sparce matrix oracle
use std::fmt::Debug;
use std::ops::{Neg, AddAssign, Mul};
use std::cmp::PartialEq;
use ndarray::Array;
use ndarray::array;
use ndarray_npy::write_npy;
use std::io; 
use ordered_float::OrderedFloat;
mod tri_loss_gurobi;
use tri_loss_gurobi::*;
mod edge_loss_gurobi;                
use edge_loss_gurobi::*;
type Coefficient = Ratio<i16>;


fn main() {    
    let mut f = BufReader::new(File::open("data_text/gamma-4-dis_mat.txt").unwrap());
    let mut s = String::new();

     // for the input as Vec-of-Vec square symmetric matrix
     let arr: Vec<Vec<f64>> = f.lines()
     .map(|l| l.unwrap().split(char::is_whitespace)
          .map(|number| number.parse().unwrap())
          .collect())
     .collect();
    let dismat = tri_loss_gurobi::ordered_floats_nested(arr.clone());
    
    let mut max_threshold = f64::INFINITY;
    for row in arr.clone(){
        let row_max = row.iter().cloned().fold(0./0., f64::max);
        if row_max < max_threshold{
            max_threshold = row_max;
        }
    }

    // set the max dimension to compute persistent homology
    let mut dim_input = String::new();
    println!("Set the dimension of the persistent homology:");
    let mut dim = 1;
    std::io::stdin().read_line(&mut dim_input);
    dim = dim_input.trim().parse::<usize>().unwrap();

    // set the maximum dissimilarity threshold
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
    
    // obtain a list of (birth_edge, death_triangle) pairs for the nonzero bars 
    let simplex_bar = simplex_barcode( &factored_complex, dim );
    
    // Get input for which method to use
    let mut is_tri_input = String::new();
    println!("Are you going to use the edge loss method? (y/n)");
    let mut is_tri= false;
    std::io::stdin().read_line(&mut is_tri_input);
    if is_tri_input.trim()=="n"{
        is_tri = true;
    }

    // Get input for whether the coefficient is integer
    let mut is_int_input = String::new();
    println!("Are the coefficients integer? (y/n)");
    let mut is_int = true;
    std::io::stdin().read_line(&mut is_int_input);
    if is_int_input.trim()=="n"{
        is_int = false;
    }
    
    // Get input for the weight method
    let mut is_uniform_input = String::new();
    println!("Is the weight method uniform? (y/n)");
    let mut is_uniform = true;
    std::io::stdin().read_line(&mut is_uniform_input);
    if is_uniform_input.trim()=="n"{
        is_uniform = false;
    }

    // Get input for the folder name 
    let mut folder_name_input = String::new();
    println!("Give the name of the folder that will store all the result data:");
    std::io::stdin().read_line(&mut folder_name_input);
    folder_name_input = folder_name_input.trim().to_string();
    
    // Set the index of the cycle that needs to be optimized
    let mut indices_input = String::new();
    println!("Which cycle(s) do you want to optimize?
Input the index of the cycle as a number or a list of numbers separated by whitespace. 
Enter \"all\" if you'd like to optimize all the cycles:");
    let mut indices = Vec::new();
    let mut is_all = false;
    std::io::stdin().read_line(&mut indices_input);
    if indices_input.trim()=="all" {
        is_all = true;
    }
    else{
        indices = indices_input.split_whitespace()
        .map(|s| s.parse().expect("parse error"))
        .collect();   
    }

    
    for j in 0..simplex_bar.len(){
        if (!is_all) && (!indices.contains(&j)){
            continue;
        }

        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        // println!("birth: {:?} death: {:?}",birth,death);
        // Write solution to npy
        
        let mut solution_hash_map = HashMap::new();
        if is_tri{
            let mut prev_obj_value = INFINITY;
            if is_uniform{
                let result_pos = tri_opt(&factored_complex, birth,death, dim,is_int,true, |x| 1.0,prev_obj_value);
                solution_hash_map = result_pos.0;
                prev_obj_value = result_pos.1;
                let result_neg = tri_opt(&factored_complex, birth,death, dim,is_int,false, |x| 1.0,prev_obj_value);
                if result_neg.0.len() > 0{
                    solution_hash_map = result_neg.0;
                }
            }
            else{
                let result_pos = tri_opt(&factored_complex, birth,death, dim,is_int,true, |x| getArea(&x, &dismat),prev_obj_value);
                solution_hash_map = result_pos.0;
                prev_obj_value = result_pos.1;
                let result_neg = tri_opt(&factored_complex, birth,death, dim,is_int,false, |x| getArea(&x, &dismat),prev_obj_value);
                if result_neg.0.len() > 0{
                    solution_hash_map = result_neg.0;
                }
            }
        }else{
            if is_uniform{
                solution_hash_map = edge_opt(&factored_complex, birth,death, dim,is_int, |x| 1.0);
            }
            else{
                solution_hash_map = edge_opt(&factored_complex, birth,death, dim,is_int, |x| getLength(&x, &dismat));
            }
        }

        
        

        let mut vertices_sol_vec = Vec::new();
        let mut coeff_sol_vec = Vec::new();
        // println!("weight {:?}",solution_hash_tri);

        for (print_key, print_val) in solution_hash_map.iter() {
            vertices_sol_vec.push(print_key.vertices[0]);
            vertices_sol_vec.push(print_key.vertices[1]);
            coeff_sol_vec.push(*print_val);
        }

        let vertices_sol_arr = Array::from_vec(vertices_sol_vec);
        let coeff_sol_arr = Array::from_vec(coeff_sol_vec); 
        
        // Create a folder to hold the results
        let folder_path = format!("{}{}{}{}", "/Users/26389/github/optimal_reps/",folder_name_input,"_", j);
        fs::create_dir(folder_path.clone());
        
        if is_tri{
            if is_uniform{
                write_npy(folder_path.clone() + "/uniform_tri_loss_vertices.npy", &vertices_sol_arr);
                write_npy(folder_path.clone() + "/uniform_tri_loss_coeffs.npy", &coeff_sol_arr);
            } 
            else{
                write_npy(folder_path.clone() + "/weighted_tri_loss_vertices.npy", &vertices_sol_arr);
                write_npy(folder_path.clone() + "/weighted_tri_loss_coeffs.npy", &coeff_sol_arr);
            }
        }
        else{
            if is_uniform{
                write_npy(folder_path.clone() + "/uniform_edge_loss_vertices.npy", &vertices_sol_arr);
                write_npy(folder_path.clone() + "/uniform_edge_loss_coeffs.npy", &coeff_sol_arr);
            } 
            else{
                write_npy(folder_path.clone() + "/weighted_edge_loss_vertices.npy", &vertices_sol_arr);
                write_npy(folder_path.clone() + "/weighted_edge_loss_coeffs.npy", &coeff_sol_arr);
            }
        }
        


        // // Write original basis to npy
    //     let x_orig = factored_complex.get_matched_basis_vector(1, &birth);
    //     println!("orig {:?}",x_orig);
    //     let mut vertices_orig_vec = Vec::new();
    //     let mut coeff_orig_vec: std::vec::Vec::<f64> = Vec::new();

    //     for (print_key, print_val) in x_orig.iter() {
    //         vertices_orig_vec.push(print_key.vertices[0]);
    //         vertices_orig_vec.push(print_key.vertices[1]);
    //         coeff_orig_vec.push((print_val.numer()/print_val.denom()).into());
    //     }
    //     let vertices_orig_arr = Array::from_vec(vertices_orig_vec);
    //     let coeff_orig_arr = Array::from_vec(coeff_orig_vec);

    //     write_npy("npy_files_dis_mat_cycle_2/orig_vertices.npy", &vertices_orig_arr);
    //     write_npy("npy_files_dis_mat_cycle_2/orig_coeffs.npy", &coeff_orig_arr);
    }
}