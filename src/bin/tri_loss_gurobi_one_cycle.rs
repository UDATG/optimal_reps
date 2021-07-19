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
use exhact::solver::multiply_hash_smoracle; // multiply a hashmap by a sparce matrix oracle
use std::fmt::Debug;
use std::ops::{Neg, AddAssign, Mul};
use std::cmp::PartialEq;
use ndarray::Array;
use ndarray::array;
use ndarray_npy::write_npy;
use ordered_float::OrderedFloat;

type Coefficient = Ratio<i16>;

fn ordered_floats( v : Vec<f64> ) -> Vec< OrderedFloat<f64> > {
    let u : Vec<_> = v.into_iter().map(OrderedFloat).collect(); 
    return u
}

fn ordered_floats_nested(v: Vec<Vec<f64>>) -> Vec< Vec< OrderedFloat<f64> > > {
    return v.into_iter().map( ordered_floats ).collect();
}

fn rational_to_float(r: Ratio<i16>)->f64{
    return r.numer().clone() as f64/r.denom().clone() as f64;
}


fn getArea( simp: &Simplex<OrderedFloat<f64>>, dismat: &Vec<Vec<OrderedFloat<f64>>> ) -> f64 {

    
   let a = f64::from(dismat[simp.vertices[0] as usize][simp.vertices[1] as usize]);
   let b = f64::from(dismat[simp.vertices[0] as usize][simp.vertices[2] as usize]);
   let c = f64::from(dismat[simp.vertices[1] as usize][simp.vertices[2] as usize]); // bug?
   let s = (a + b + c)/2.;
   let t = s*(s-a)*(s-b)*(s-c);
   if t<0.{
       println!("Triangle inequalty is violated. This message is generated in function getArea.");
   }
   return t.sqrt();
 }

fn tri_opt<'a, MatrixIndexKey, Filtration, OriginalChx, Matrix, WeightFunction>
(   factored_complex: &FactoredComplexBlockCsm<'a, MatrixIndexKey, Coefficient, Filtration, OriginalChx>,
    birth: &MatrixIndexKey, 
    death: &MatrixIndexKey,
    dim: usize, 
    is_int: bool, 
    is_pos: bool,
    weight: WeightFunction )
    
    -> HashMap<MatrixIndexKey, f64> 
    where   OriginalChx: ChainComplex<MatrixIndexKey, Coefficient, Filtration, Matrix=Matrix>,
            MatrixIndexKey: PartialEq+ Eq + Clone + Hash + std::cmp::PartialOrd + Ord + std::fmt::Debug,
            Matrix: SmOracle<MatrixIndexKey, MatrixIndexKey, Coefficient>,
            Filtration: PartialOrd + Clone,
            WeightFunction: Fn( &MatrixIndexKey ) -> f64{
    
        let chx = &factored_complex.original_complex;
        let i = dim;
        let env = Env::new("logfile.log").unwrap();

        let mut model = env.new_model("model_pos").unwrap();// D_{n+1}[\sigma, \hat{F}_{n+1}]v>0
    
        // a list of tuples (birth simplex, death simplex)
        // loop over Sn+1
        //replace inequality with strict inequality
        let mut good_triangles: Vec<_> = chx.keys_unordered_itr(i + 1).filter(|s| chx.key_2_filtration(&s) >= chx.key_2_filtration(&birth) && &s <= &death ).collect();
        
        
        // count the size of good triangles
        // println!("good triangle:{:?}",good_triangles);

        let size = good_triangles.len();
        
        
         // build oracle for the entire boundary matrix
        let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row, exhact::chx::ChxTransformKind::Boundary);


        // Create F_n
        // sigma'>sigma in the linear order
        let mut good_edges: Vec<_> = chx.keys_unordered_itr(i).filter(|s| chx.key_2_filtration(&s) <= chx.key_2_filtration(&death) && &s > &birth).collect();
        
        let mut death_iter = D.min_itr(&death);
        for death_edge in death_iter{
            println!("death edge:{:?}", death_edge);
            println!("death edge in good edge :{:?}", good_edges.contains(&death_edge.0));
        }

        // println!("good edges:{:?}", good_edges);
        // Set program type
        let mut program_type = Integer;
        if (is_int){
            program_type = Integer;
        }
        else{
            program_type = Continuous;
        }


        // create hashmaps to store triangles to indices 
        let mut triangle_2_index: HashMap<MatrixIndexKey, usize> = HashMap::new();       
        let mut index_2_triangle: HashMap<usize, MatrixIndexKey> = HashMap::new();       
        // initialize indices to be 0   
        let mut maj_index:usize = 0;
        
        for triangle in good_triangles { // for each column
            if !triangle_2_index.contains_key(&triangle) {
                triangle_2_index.insert(triangle.clone(), maj_index.clone());
                index_2_triangle.insert(maj_index.clone(), triangle.clone());
                maj_index = maj_index + 1;
            }
        }
        

        // set up variables v = (v+)-(v-)
        // optimize for variables [v+ v-]
        // set up the vector: v+
        let mut v_pos = Vec::new();
        for i in 0..size{
            let mut name = format!("{}{}", "v_pos", i);
            let mut str_name = &name[..];
            // set up for variables in model_pos and model_neg
            // variable range: 0-INFINITY
            v_pos.push(model.add_var(str_name, program_type, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
        }
        // set the vector: v-
        let mut v_neg = Vec::new();
        for i in 0..size{
            let mut name = format!("{}{}", "v_neg", i);
            let mut str_name = &name[..];
            // set up for variables in model_pos and model_neg
            // variable range: 0-INFINITY
            v_neg.push(model.add_var(str_name, program_type,1.0, 0.0, INFINITY, &[], &[]).unwrap());
        }
        
        // Set objective function
        let mut obj_expression: LinExpr = LinExpr::new();
        for i in 0..size{    
            obj_expression = obj_expression.add_term(weight(index_2_triangle.get(&i).unwrap()), v_pos[i].clone());
        }
        for i in 0..size{    
            obj_expression = obj_expression.add_term(weight(index_2_triangle.get(&i).unwrap()), v_neg[i].clone());
        }
      

         // integrate all of the constraints into the model_pos.
         model.update().unwrap();
         






        // Set up constraints 

        //Set up the first kind of constraint
        // D_{n+1}[sigma,\hat{F}_{n+1}] v != 0 
        let row = D.maj_itr(&birth);
        let mut constraint1 = LinExpr::new();
        for item in row{
            // Check whether the item is in the set of good triangels, which corresponds to F_{n+1}_hat in the paper
            // This differs from F_{n+1} by addong the death simplex

            // if (chx.key_2_filtration(&item.0) >= chx.key_2_filtration(&birth) && &item.0 < &death){
            if (&item.0 <= &death){
                let mut index: usize = triangle_2_index.get(&item.0).unwrap().clone();
                // set coefficients
                constraint1 = constraint1.add_term((item.1.numer()/item.1.denom()).into(), v_pos[index].clone());
                constraint1 = constraint1.add_term((-item.1.numer()/item.1.denom()).into(), v_neg[index].clone());
            }
        }
        // set constraintSense to be either greater or less than 0
        if is_pos{
            model.add_constr("constraint1", constraint1.clone(), Greater, f64::EPSILON);
        }
        else{
            model.add_constr("constraint1", constraint1.clone(), Less, -f64::EPSILON);
        }
        
        
        model.update().unwrap();
        

        //Set up the second kind of constraint
        // D_{n+1}[F_n,\hat{F}_{n+1}] v == 0 
        for edge in good_edges{
            //println!("{:?} good_edges", edge.clone());
            
            let mut row_ctr2 = D.maj_itr(&edge);
            let mut constraint2 = LinExpr::new();
            for item in row_ctr2{
                // Check whether the item is in the set of good triangels, which corresponds to F_{n+1}_hat in the paper
                // This differs from F_{n+1} by addong the death simplex
                
                // if (chx.key_2_filtration(&item.0) >= chx.key_2_filtration(&birth) && &item.0 < &death){
                if (&item.0 <= &death){
                    let mut index: usize = triangle_2_index.get(&item.0).unwrap().clone();
                    constraint2 = constraint2.add_term((item.1.numer()/item.1.denom()).into(), v_pos[index].clone());
                    constraint2 = constraint2.add_term((-item.1.numer()/item.1.denom()).into(), v_neg[index].clone());
                }
            }
            // if (edge. == 1 && edge.vertices[1] == 67){
            //     println!("wrong constraint {:?}",constraint2);
            // }
            model.add_constr("constraint2", constraint2.clone(), Equal, 0.0);
            model.update().unwrap();
        }
        
        // Set up the third kind of constraint
        // v (death) = 1
        let mut constraint3 = LinExpr::new();
        let mut index_ctr3: usize = triangle_2_index.get(&death).unwrap().clone();
        constraint3 = constraint3.add_term(1.0,v_pos[index_ctr3].clone());
        constraint3 = constraint3.add_term(-1.0,v_neg[index_ctr3].clone());
        model.add_constr("constraint3", constraint3.clone() , Equal, 1.0);
        
        

        model.update().unwrap();
        

        // add objective function to model_pos
        model.set_objective(obj_expression.clone(),Minimize).unwrap();
        model.update().unwrap();
        model.write("logfile_tri_loss_new.lp").unwrap();
        model.optimize().unwrap();
        
        let v_neg_val = model.get_values(attr::X, &v_neg).unwrap(); // Get the result for x+ and x-
        let v_pos_val = model.get_values(attr::X, &v_pos).unwrap();

    //     let mut ans: HashMap<MatrixIndexKey, f64> = HashMap::new();
    //     for i in 0..size{
    //     if v_pos_val[i]- v_neg_val[i]!=0.{
    //         let tri_ans = index_2_triangle.get(&i).unwrap().clone();
    //         let tri_val = v_pos_val[i]- v_neg_val[i];
    //         ans.insert(tri_ans, tri_val);
    //         }
    //     }

    // return ans;
    let mut solution_hash_edge : HashMap<MatrixIndexKey, f64> = HashMap::new();  

        // if model_pos.get(attr::ObjVal).unwrap()<=model_neg.get(attr::ObjVal).unwrap(){ // if model_pos gives the smaller objective value
            let v_neg_val = model.get_values(attr::X, &v_neg).unwrap();
            let v_pos_val = model.get_values(attr::X, &v_pos).unwrap();
            // hashmap (triangle, coefficients) for triangles present in the 2-chain
            let mut solution_hash_triangle: HashMap<MatrixIndexKey, f64> = HashMap::new();  
            for i in 0..size{
                let sol_value = v_pos_val[i]- v_neg_val[i];
                if sol_value!=0.{
                    solution_hash_triangle.insert(index_2_triangle.get(&i).unwrap().clone(), sol_value);
                }
            }
            println!("tri ans length : {:?}",solution_hash_triangle.len());
            let top_tri = solution_hash_triangle.keys().max().clone();
            println!("top tri:{:?}",top_tri);
            // build solution_hash_edge
            for (tri_key, tri_val) in solution_hash_triangle {
                for ( edge_key, edge_val) in D.min_itr(& tri_key ) {
                    if !solution_hash_edge.contains_key(&edge_key){
                        // add key 
                        // val = edge_val*tri_val
                        solution_hash_edge.insert(edge_key, rational_to_float(edge_val)*tri_val);
                    }
                    else{
                        // update val = current_val+edge_val*tri_val
                        let mut current_val = solution_hash_edge.get_mut(&edge_key).unwrap();
                        *current_val =  *current_val+(rational_to_float(edge_val))*tri_val;
                    }
                }
            }

            let support : Vec<_> = solution_hash_edge.keys().cloned().collect();
            for edge_key in support{
                if solution_hash_edge.get(&edge_key).unwrap()==&0.{
                    solution_hash_edge.remove(&edge_key);
                }
            }
        // }
        // else{ // else if model_neg gives the smaller objective value
        //     let v_neg_val = model_neg.get_values(attr::X, &v_neg).unwrap();
        //     let v_pos_val = model_neg.get_values(attr::X, &v_pos).unwrap();
        //     // hashmap (triangle, coefficients) for triangles present in the 2-chain
        //     let mut solution_hash_triangle: HashMap<MatrixIndexKey, f64> = HashMap::new();  
        //     for i in 0..size{
        //         let sol_value = v_pos_val[i]- v_neg_val[i];
        //         if sol_value!=0.{
        //             solution_hash_triangle.insert(index_2_triangle.get(&i).unwrap().clone(), sol_value);
        //         }
        //     }
        //     println!("tri ans length: {:?}",solution_hash_triangle.len());
        //     // build solution_hash_edge
        //     for (tri_key, tri_val) in solution_hash_triangle {
        //         for ( edge_key, edge_val) in D.min_itr(& tri_key ) {
        //             if !solution_hash_edge.contains_key(&edge_key){
        //                 // add key 
        //                 // val = edge_val*tri_val
        //                 solution_hash_edge.insert(edge_key, rational_to_float(edge_val)*tri_val);
        //             }
        //             else{
        //                 // update val = current_val+edge_val*tri_val
        //                 let mut current_val = solution_hash_edge.get_mut(&edge_key).unwrap();
        //                 *current_val =  *current_val+(rational_to_float(edge_val))*tri_val;
        //             }
        //         }
        //     }

        //     let support : Vec<_> = solution_hash_edge.keys().cloned().collect();
        //     for edge_key in support{
        //         if solution_hash_edge.get(&edge_key).unwrap()==&0.{
        //             solution_hash_edge.remove(&edge_key);
        //         }
        //     }

        // }
        return solution_hash_edge;


        
}

fn main() {    
    // read distance matrix
    let mut f = BufReader::new(File::open("data_text\\dist_mat.txt").unwrap());
    let mut s = String::new();

     // for the input as Vec-of-Vec square symmetric matrix
     let arr: Vec<Vec<f64>> = f.lines()
     .map(|l| l.unwrap().split(char::is_whitespace)
          .map(|number| number.parse().unwrap())
          .collect())
     .collect();
    let dismat = ordered_floats_nested(arr);
    
    
    // set the max dimension to compute persistent homology
    let dim = 1;

    // set the maximum dissimilarity threshold
    const maxdis: OrderedFloat<f64> = OrderedFloat(4.);

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
    let simplex_bar = simplex_barcode( &factored_complex, 1 );
    let cycle_index = 2;

    for j in cycle_index..cycle_index+1{
        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        // println!("birth: {:?} death: {:?}",birth,death);
        // Write solution to npy
        
        let solution_hash_tri = tri_opt(&factored_complex, birth,death, 1,true,false, |x| getArea(&x, &dismat));



        let mut vertices_sol_vec = Vec::new();
        let mut coeff_sol_vec = Vec::new();
        // let mut sum:f64 = 0.0;
        for (print_key, print_val) in solution_hash_tri.iter() {
            vertices_sol_vec.push(print_key.vertices[0]);
            vertices_sol_vec.push(print_key.vertices[1]);
            coeff_sol_vec.push(*print_val);

            // sum += chx.key_2_filtration(&print_key).into_inner();
        }
        // println!("sum : {:?}",sum);

        let vertices_sol_arr = Array::from_vec(vertices_sol_vec);
        let coeff_sol_arr = Array::from_vec(coeff_sol_vec);
        let folder_name = format!("{}{}", "npy_files_dis_mat_cycle_tri", j);  
        write_npy(folder_name.clone() + "/weighted_edge_answer_vertices.npy", &vertices_sol_arr);
        write_npy(folder_name.clone() + "/weighted_edge_answer_coeffs.npy", &coeff_sol_arr);
        println!("{}",folder_name);


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