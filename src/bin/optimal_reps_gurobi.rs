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

// must add the following line to dependencies under Cargo.toml:
//      ordered-float = "2.0"

use ordered_float::OrderedFloat;
type Coefficient = Ratio<i16>;

fn ordered_floats( v : Vec<f64> ) -> Vec< OrderedFloat<f64> > {
    let u : Vec<_> = v.into_iter().map(OrderedFloat).collect(); 
    return u
}

fn ordered_floats_nested(v: Vec<Vec<f64>>) -> Vec< Vec< OrderedFloat<f64> > > {
    return v.into_iter().map( ordered_floats ).collect();
}

pub enum SimplexWeights {
    Uniform,
    Volume,
}

// pub enum ProgramType {
//     LP,
//     MIP,
// }

fn tri_opt<'a, MatrixIndexKey, Filtration, OriginalChx, Matrix>
(
    is_pos: bool, // optimize over the positive domain
    is_int: bool,
    dim: usize,
    weight: SimplexWeights, 
    factored_complex: &FactoredComplexBlockCsm<'a, MatrixIndexKey, Coefficient, Filtration, OriginalChx>,
    birth: &MatrixIndexKey,
    death: &MatrixIndexKey)//-> sprs::CsVecBase<std::vec::Vec<usize>, std::vec::Vec<f64>, f64> // Vec<f64>
    where   OriginalChx: ChainComplex<MatrixIndexKey, Coefficient, Filtration, Matrix=Matrix>,
            MatrixIndexKey: PartialEq+ Eq + Clone + Hash + std::cmp::PartialOrd + Ord + std::fmt::Debug,
            Matrix: SmOracle<MatrixIndexKey, MatrixIndexKey, Coefficient>,
            Filtration: PartialOrd + Clone,
    {
        let chx = &factored_complex.original_complex;
        let i = dim;
        let env = Env::new("logfile.log").unwrap();
        let mut model = env.new_model("model1").unwrap();
       
        println!("{:?}", birth);
        println!("{:?}", death);
        // a list of tuples (birth simplex, death simplex)
        // loop over Sn+1
        let good_triangles = chx.keys_unordered_itr(i + 1).filter(|s| s <= &death && s >= &birth );
        let size = good_triangles.count();
        // loop over Sn
        let good_edges = chx.keys_unordered_itr(i).filter(|s| s <= &death && s >= &birth);
        let obj_coef;
        match weight {
            SimplexWeights::Uniform => {
                obj_coef = vec![1.; 2 * size]; // c^T // 1 vector with length |Fn|
            }
            SimplexWeights::Volume => {
                obj_coef = vec![1.; 2 * size]; // c^T // 1 vector with length |Fn|
            }
        }
        
        let mut program_type = Integer;
        if (is_int){
            program_type = Integer;
        }
        else{
            program_type = Continuous;
        }
        

        // initialize the vector: v+
        let mut v_pos = Vec::new();
        
        for i in 0..size{

            let mut name = format!("{}{}", "x_pos", i);

            let mut str_name = &name[..];
            v_pos.push(model.add_var(str_name, program_type, 0.0, -INFINITY, INFINITY, &[], &[]).unwrap());
        }

        // initialize the vector: v-
        let mut v_neg = Vec::new();
        
        for i in 0..size{

            let mut name = format!("{}{}", "x_neg", i);

            let mut str_name = &name[..];
            v_neg.push(model.add_var(str_name, program_type, 0.0, -INFINITY, INFINITY, &[], &[]).unwrap());
        }
        
        // Set objective function
        let mut obj_expression: LinExpr = LinExpr::new();
        

        let coeffs_vec= vec![1.0; size]; 
        let coeffs: &[f64] = &coeffs_vec[..];

        let a:&[f64] =&[1.0,11.0,1.0];
        // obj_expression.add_terms(coeffs: a, vars: &v_pos[..]);

        // for i in 0..size{    
            
        //     obj_expression.add_term(1.0, v_pos[i].clone());

        // }



//             // build oracle for the entire boundary matrix
        let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row, exhact::chx::ChxTransformKind::Boundary);
        let maj_fields = D.min_itr(&death);
        for item in maj_fields{
            println!("{:?}",item);
            println!("{}", &item.0 >= &birth && &item.0 <= &death);
        }
//             // create hashmaps to store keys to indices 
//             let mut maj_2_index: HashMap<MatrixIndexKey, usize> = HashMap::new();       
//             let mut min_2_index: HashMap<MatrixIndexKey, usize> = HashMap::new();   
//             let mut index_2_maj: HashMap<usize, MatrixIndexKey> = HashMap::new();       
//             let mut index_2_min: HashMap<usize, MatrixIndexKey> = HashMap::new();  
//             // initialize indices to be 0   
//             let mut maj_index:usize = 0;
//             let mut min_index:usize = 0;
//             // create sparse matrix
//             let mut ind_ptr = Vec::new();
//             ind_ptr.push(0);
//             let mut col_ind = Vec::new();
//             let mut nz_val : Vec<f64> = Vec::new();
//             let mut counter = 1;
//             for edge in Fn { // for each row 
//                 // println!("{}", counter);
//                 let row = m.add_row();
//                 if &edge == birth {
//                     if is_pos{
//                         m.set_row_lower(row, f64::EPSILON);
//                     }
//                     else{
//                         m.set_row_upper(row, -f64::EPSILON);
//                     }
//                 }
//                 else{
//                     println!("0000");
//                     println!("{:?}", edge);
//                     m.set_row_upper(row, 0.0);
//                     m.set_row_lower(row, 0.0);
//                 }
//                 if !maj_2_index.contains_key(&edge) {
//                     maj_2_index.insert(edge.clone(), maj_index.clone());
//                     index_2_maj.insert(maj_index.clone(), edge.clone());
//                     maj_index = maj_index + 1;
//                 }
//                 let minor_fields = D.maj_itr(&edge);
//                 for minor_field in minor_fields{
//                     // column index (S_{n+1})
//                     let tri = minor_field.0;
//                     // entry value
//                     let data = minor_field.1;
                    
//                     // if &tri == death{
//                     //     println!("tau");
//                     //     println!("{:?}", min_2_index.contains_key(&tri));
//                     //     println!("the edge that's contained in tau");
//                     //     println!("{:?}", edge);
//                     // }
//                     if &tri <= &death && &tri >= &birth{                        
//                         if !min_2_index.contains_key(&tri) {
//                             min_2_index.insert(tri.clone(), min_index.clone());
//                             index_2_min.insert(min_index.clone(), tri.clone());
//                             min_index = min_index + 1;
//                         }
//                         col_ind.push(*min_2_index.get(&tri).unwrap());
//                         nz_val.push((data.numer()/data.denom()).into());
//                         counter = counter+1;
//                         m.set_weight(row, cols[*min_2_index.get(&tri).unwrap()], (data.numer()/data.denom()).into());
//                         m.set_weight(row, cols[*min_2_index.get(&tri).unwrap() + size], (-data.numer()/data.denom()).into());
//                     }
                    
//                 }
//                 ind_ptr.push(counter);
//             }
//             // let a = CsMat::new((3, 3),
//             //            vec![0, 2, 4, 5],
//             //            vec![0, 1, 0, 2, 2],
//             //            vec![1., 2., 3., 4., 5.]);
                       
//             // println!("{}", ind_ptr.clone().nnz());
//             // println!("{}", nz_val.clone().len());
//             // let csr = CsMat::new((ind_ptr.clone().len() - 1, size), ind_ptr, col_ind, nz_val);
//             // Set objective sense.
//             m.set_obj_sense(Sense::Minimize);
//             m.set_col_upper(cols[*min_2_index.get(&death).unwrap()], 1.0);
//             m.set_col_lower(cols[*min_2_index.get(&death).unwrap()], 1.0);
//             m.set_col_upper(cols[*min_2_index.get(&death).unwrap() + size], 0.0);
//             m.set_col_lower(cols[*min_2_index.get(&death).unwrap() + size], 0.0);
//             // Solve the problem. Returns the solution
//             let sol = m.solve();
            
//             let mut ind = Vec::new();
//             let mut val = Vec::new();
//             for i in 0..size{
//                 if sol.col(cols[i]) - sol.col(cols[i + size]) != 0. {
//                     ind.push(i);
//                     val.push(sol.col(cols[i]) - sol.col(cols[i + size]));
//                     println!("{:?}", index_2_min.get(&i));
//                 }
//             }
//             let mut v = CsVec::new(size, ind, val);
//             // let x = &csr * &v;
//             println!("{:?}", v);
//         return v;
}

fn main() {    
    let mut f = BufReader::new(File::open("senate104_edge_list.txt_0.68902_distmat.txt").unwrap());
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
    const maxdis: OrderedFloat<f64> = OrderedFloat(1.);

    // create a "ring object" representing the field of rational numbers
    let ringmetadata = exhact::matrix::RingMetadata{
    ringspec: RingSpec::Rational,
    identity_additive: Ratio::new(0, 1),
    identity_multiplicative: Ratio::new(1, 1)
    };

    println!("here2");
    // build and factor the filtered chain complex
    let chx = exhact::clique::CliqueComplex {
        // the distance/dissimilarity matrix
        dissimilarity_matrix: dismat, 
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
    println!("here3");
    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);
    
    println!("here4");
    // obtain a list of (birth_edge, death_triangle) pairs for the nonzero bars 
    let simplex_bar = simplex_barcode( &factored_complex, 1 );
    

    for j in 1..2{//simplex_bar.len(){
        
        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        let v = tri_opt(true,true, 1, SimplexWeights::Uniform, &factored_complex, birth, death);
        //sols.push(v);
    }

}
