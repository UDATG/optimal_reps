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

//skip alg 1, focus on program 14

// For set Q, take pivot column with birth time smaller than x_orig.

type Coefficient = Ratio<i16>;
use ordered_float::OrderedFloat;

fn ordered_floats( v : Vec<f64> ) -> Vec< OrderedFloat<f64> > {
    let u : Vec<_> = v.into_iter().map(OrderedFloat).collect(); 
    return u
}

fn ordered_floats_nested(v: Vec<Vec<f64>>) -> Vec< Vec< OrderedFloat<f64> > > {
    return v.into_iter().map( ordered_floats ).collect();
}


fn edge_opt<'a, MatrixIndexKey, Filtration, OriginalChx, Matrix>(
    factored_complex: &FactoredComplexBlockCsm<'a, MatrixIndexKey, Coefficient, Filtration, OriginalChx>,
    birth: &MatrixIndexKey, death: &MatrixIndexKey, dim: usize, is_int:bool, is_uniform:bool
) ->  HashMap<MatrixIndexKey, f64>
where   OriginalChx: ChainComplex<MatrixIndexKey, Coefficient, Filtration, Matrix=Matrix>,
            MatrixIndexKey: PartialEq+ Eq + Clone + Hash + std::cmp::PartialOrd + Ord + std::fmt::Debug,
            Matrix: SmOracle<MatrixIndexKey, MatrixIndexKey, Coefficient>,
            Filtration: PartialOrd + Clone,
{

    let chx = &factored_complex.original_complex;
    let i = dim;
    let env = Env::new("logfile.log").unwrap();
    let mut model = env.new_model("model1").unwrap();

    

    let good_edges: Vec<_> = chx.keys_unordered_itr(i).filter(|s| s <= &birth).collect();
    println!("test");
    
    let mut good_triangles:Vec<_> = Vec::new();
    
    
    let simplex_bar = simplex_barcode( &factored_complex, 1 );
    let barcode_length = simplex_bar.len();

    // Build matrix A and good triangles at the same time
    let mut A = Vec::new();
    for bar_index in 0..barcode_length{
        let  barcode_birth = &simplex_bar[bar_index].0;
        let  barcode_death = &simplex_bar[bar_index].1;

        // Build good triangle
        if barcode_death <= birth{
            good_triangles.push(barcode_death);
        }
        
        // Build matrix A
        if barcode_birth <= birth && barcode_death <= death{
            A.push(factored_complex.get_matched_basis_vector(1, &barcode_birth));
        }
    }
    let edge_size = good_edges.len();
    let triangle_size = good_triangles.len();
    let column_size_of_a = A.len(); 

    let mut weight = Vec::new();
    //Set weight type
    if is_uniform{
        for i in 0..edge_size{
            weight.push(1.0);
        }
    }
    else{
        for edge in good_edges.iter(){
            // let edge_length = chx.key_2_filtration(&edge).into_inner();

            // weight.push(edge_length);
        }
    }
    

    // Set program type
    let  mut program_type = Integer;
    if is_int{
        program_type = Integer;
    }
    else{
        program_type = Continuous;
    }

    // initialize the vector: x+
    let mut x_pos = Vec::new();
        
    for i in 0..edge_size{

        let  name = format!("{}{}", "x_pos", i);

        let  str_name = &name[..];
        x_pos.push(model.add_var(str_name, program_type, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: x-
    let mut x_neg = Vec::new();
    
    for i in 0..edge_size{

        let  name = format!("{}{}", "x_neg", i);

        let  str_name = &name[..];
        x_neg.push(model.add_var(str_name, program_type, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: p 
    let mut p = Vec::new();
    
    for i in 0..column_size_of_a{

        let  name = format!("{}{}", "p", i);

        let  str_name = &name[..];
        p.push(model.add_var(str_name, program_type, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: q
    let mut q = Vec::new();
    
    for i in 0..triangle_size{

        let  name = format!("{}{}", "q", i);

        let  str_name = &name[..];
        q.push(model.add_var(str_name, program_type, 1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // Set objective function
    let mut obj_expression: LinExpr = LinExpr::new();
        
    for i in 0..edge_size{    
        
        obj_expression = obj_expression.add_term(weight[i], x_pos[i].clone());

    }

    for i in 0..edge_size{    
        
        obj_expression = obj_expression.add_term(weight[i], x_neg[i].clone());

    }
    


    // The cycle vector: x_orig
    let x_orig = factored_complex.get_matched_basis_vector(1, &birth);
    
    // Build HashMap to record the index of the edges
    let mut edge_2_index: HashMap<MatrixIndexKey, usize> = HashMap::new();       
    let mut index_2_edge: HashMap<usize, MatrixIndexKey> = HashMap::new();
    let mut edge_index:usize = 0;
    for edge in good_edges.iter(){
        if !edge_2_index.contains_key(edge){
            edge_2_index.insert(edge.clone(), edge_index.clone());
            index_2_edge.insert(edge_index.clone(), edge.clone());
            edge_index = edge_index + 1;
        }
    }

    // Build Hashmap to record the index of the triangles
    let mut triangle_2_index: HashMap<&MatrixIndexKey, usize> = HashMap::new();       
    let mut index_2_triangle: HashMap<usize, &MatrixIndexKey> = HashMap::new();       
    // initialize indices to be 0   
    let mut tri_index:usize = 0;
    
    // How do we build good_triangles?


    for triangle in good_triangles.iter() { // for each column
        if !triangle_2_index.contains_key(triangle) {
            triangle_2_index.insert(triangle.clone(), tri_index.clone());
            index_2_triangle.insert(tri_index.clone(), triangle.clone());
            tri_index = tri_index + 1;
        }
    }

    let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row, exhact::chx::ChxTransformKind::Boundary);

    // Set constraint
    for edge in good_edges.iter(){
        let mut row = D.maj_itr(edge);
        let mut constraint = LinExpr::new();
        // Add x+ and x-
        let index = edge_2_index.get(&edge).unwrap().clone();
        constraint = constraint.add_term(1.0, x_pos[index].clone());
        constraint = constraint.add_term(-1.0, x_neg[index].clone());

        // Add a constant which represents x_orig
        if x_orig.contains_key(&edge){
            let x_orig_val = x_orig.get(&edge).unwrap().clone(); 
            constraint = constraint.add_constant((-x_orig_val.numer()/x_orig_val.denom()).into());
        } 
        // Add a term related to vector q and the boundary matrix
        for item in row{
            if good_triangles.contains(&&item.0){
                let mut index: usize = triangle_2_index.get(&item.0).unwrap().clone();
                constraint = constraint.add_term((-item.1.numer()/item.1.denom()).into(), q[index].clone());

            }
        }

        // Add a term related to vector p and the matrix A
        for i in 0..column_size_of_a{
            let column = &A[i];
            if (column.contains_key(&edge)){
                let coeff = column.get(&edge).unwrap();
                let coeff_val = -coeff.numer()/coeff.denom();
                constraint = constraint.add_term(coeff_val.into(), p[i].clone());
            }
        }
        
        model.add_constr("constraint", constraint, Equal, 0.0);
    }

    model.update().unwrap();
    model.set_objective(obj_expression,Minimize).unwrap();
    model.write("logfile.lp").unwrap();
    model.optimize().unwrap();
    let x_neg_val = model.get_values(attr::X, &x_neg).unwrap();
    let x_pos_val = model.get_values(attr::X, &x_pos).unwrap();

    
    // let mut ind = Vec::new();
    // let mut val = Vec::new();
    let mut ans: HashMap<MatrixIndexKey, f64> = HashMap::new();
    for i in 0..edge_size{
        if x_pos_val[i]- x_neg_val[i]!=0.{
            let edge_ans = index_2_edge.get(&i).unwrap().clone();
            let edge_val = x_pos_val[i]- x_neg_val[i];
            ans.insert(edge_ans, edge_val);
        }
    }

    // let mut v = CsVec::new(edge_size, ind, val.clone());

    return ans;

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
    
    let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);
    
    
    // obtain a list of (birth_edge, death_triangle) pairs for the nonzero bars 
    let simplex_bar = simplex_barcode( &factored_complex, 1 );
    
    for j in 0..1{
        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        factored_complex.get_matched_basis_vector(1, &birth);
        let weight = &factored_complex.original_complex.key_2_filtration(&birth);
        
        let v = edge_opt(&factored_complex, birth,death, 1,true,false);
        println!("{:?}", v);
        
    }
}

