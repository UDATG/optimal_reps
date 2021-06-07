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
        
        // count the size of good triangles
        let good_triangles_copy = chx.keys_unordered_itr(i + 1).filter(|s| s <= &death && s >= &birth );
        let size = good_triangles_copy.count();

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
        // Set program type
        let mut program_type = Integer;
        if (is_int){
            program_type = Integer;
        }
        else{
            program_type = Continuous;
        }

        // Set constraint sense
        let mut constraintSense = Greater;
        if (is_pos){
            constraintSense = Greater;
        }
        else{
            constraintSense = Less;
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
        
        for i in 0..size{    
            
            obj_expression = obj_expression.add_term(1.0, v_pos[i].clone());

        }

        for i in 0..size{    
            
            obj_expression = obj_expression.add_term(1.0, v_neg[i].clone());

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


        // build oracle for the entire boundary matrix
        let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row, exhact::chx::ChxTransformKind::Boundary);

        //Set up the first kind of constraint
        let row = D.maj_itr(&birth);

        let mut constraint1 = LinExpr::new();
        for item in row{
            if (&item.0 >= &birth && &item.0 <= &death){
                let mut index: usize = *triangle_2_index.get(&item.0).unwrap();
                constraint1 = constraint1.add_term((item.1.numer()/item.1.denom()).into(), v_pos[index].clone());
                constraint1 = constraint1.add_term((-item.1.numer()/item.1.denom()).into(), v_neg[index].clone());
            }
        }

        model.add_constr("constraint1", constraint1, constraintSense, 0.0);

        //Set up the second kind of constraint
        for edge in good_edges{
            let mut row_ctr2 = D.maj_itr(&edge);
            let mut constraint2 = LinExpr::new();
            for item in row_ctr2{
                if (&item.0 >= &birth && &item.0 <= &death){
                    let mut index: usize = *triangle_2_index.get(&item.0).unwrap();
                    constraint2 = constraint2.add_term((item.1.numer()/item.1.denom()).into(), v_pos[index].clone());
                    constraint2 = constraint2.add_term((-item.1.numer()/item.1.denom()).into(), v_neg[index].clone());
                }
            }
            model.add_constr("constraint2", constraint2, Equal, 0.0);

        }

        //Set up the third kind of constraint
        // v+ (death) = 1
        let mut constraint3 = LinExpr::new();
        let mut index_ctr3: usize = *triangle_2_index.get(&death).unwrap();

        constraint3 = constraint3.add_term(1.0,v_pos[index_ctr3].clone());

        model.add_constr("constraint3", constraint3 , Equal, 1.0);

        // v- (death) = 0
        let mut constraint4 = LinExpr::new();
        let mut index_ctr4: usize = *triangle_2_index.get(&death).unwrap();

        constraint4 = constraint4.add_term(1.0,v_neg[index_ctr4].clone());

        model.add_constr("constraint4", constraint4 , Equal, 0.0);

        model.write("logfile.lp").unwrap();

        model.optimize().unwrap();
        let v_pos_val = model.get_values(attr::X, &v_pos);
        let v_neg_val = model.get_values(attr::X, &v_neg);

        println!("{:?}", v_pos_val);
        // let v = Vec::new();
        // for i in 0..size{
        //     v.push(v_pos_val[i]- v_neg_val[i]);
        // }


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
