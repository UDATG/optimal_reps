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

fn getLength( simp: &Simplex<OrderedFloat<f64>>, dismat: &Vec<Vec<OrderedFloat<f64>>> ) -> f64 {
    let a = f64::from(dismat[simp.vertices[0] as usize][simp.vertices[1] as usize]);
    return a;
  }

fn has_fifth_vertex( simp: &Simplex<OrderedFloat<f64>> ) -> f64 {
    let mut weight = 2.0;
    let constant : u16 = 4;
    for vertex in simp.vertices.iter(){
        if vertex  == &constant {
            weight = 1.0;
            break;
        }
    }
    return weight;

}


fn edge_opt<'a, MatrixIndexKey, Filtration, OriginalChx, Matrix, WeightFunction>(
    factored_complex: &FactoredComplexBlockCsm<'a, MatrixIndexKey, Coefficient, Filtration, OriginalChx>,
    birth: &MatrixIndexKey, 
    death: &MatrixIndexKey, 
    dim: usize, 
    is_int:bool, 
    weight: WeightFunction
) ->  HashMap<MatrixIndexKey, f64>
where   OriginalChx: ChainComplex<MatrixIndexKey, Coefficient, Filtration, Matrix=Matrix>,
            MatrixIndexKey: PartialEq+ Eq + Clone + Hash + std::cmp::PartialOrd + Ord + std::fmt::Debug,
            Matrix: SmOracle<MatrixIndexKey, MatrixIndexKey, Coefficient>,
            Filtration: PartialOrd + Clone ,
            WeightFunction: Fn( &MatrixIndexKey ) -> f64
{
    let i = dim;
    let chx = &factored_complex.original_complex;
    let indexing_vec = &factored_complex.dim_indexing;
    let indexing = &indexing_vec[dim+1];
    let index_2_majkey = &indexing.index_2_majkey;
    let index_2_minkey = &indexing.index_2_minkey;


    let env = Env::new("logfile.log").unwrap();
    let mut model = env.new_model("model1").unwrap();

    
    let good_edges: Vec<_> = chx.keys_unordered_itr(i).filter( |s| chx.key_2_filtration(&s) <= chx.key_2_filtration(&birth)).collect();
    let mut good_triangles:Vec<_> = Vec::new();
    // println!("good edges : {:?}", good_edges );

    // Build good triangles 
    
    for bar_index in 0..index_2_majkey.len(){
        // get the row and column index of each pivot element
        let mut index_birth = index_2_majkey[bar_index].clone(); // birth row/edge
        let mut index_death = index_2_minkey[bar_index].clone(); // death col/triangle
        
        // if the pivot has birth and death times earlier than those of the optimized cycle ..
        if chx.key_2_filtration(&index_death) <= chx.key_2_filtration(birth)
        {
            // add the column index to the set of "good triangles" when it's born before or equal to the cycle.
            // println!("index birth {:?} \n", index_birth);
            // println!("index death {:?} \n", index_death);
            good_triangles.push(index_death.clone());
            
            
        }
    }

    // Build matrix A
    let mut A = Vec::new();
    for bar_index in 0..index_2_majkey.len(){
        // get the row and column index of each pivot element
        let mut index_birth = index_2_majkey[bar_index].clone(); // birth row/edge
        let mut index_death = index_2_minkey[bar_index].clone(); // death col/triangle

        // the added class is born in interval (-inf, birth_of_optimized_class]
        // AND
        // the added class dies in interval (birth_of_optimized_class, death_of_optimized class)

        let mut condition_a =   
        chx.key_2_filtration(&index_birth)    <=  chx.key_2_filtration(&birth) &&
        chx.key_2_filtration(&index_death)    >   chx.key_2_filtration(&birth) &&
        chx.key_2_filtration(&index_death)    <   chx.key_2_filtration(&death);

        // the added class is born strictly before the birth simplex (in lexicographic order) 
        // AND
        // the added class dies when the optimized class dies

        let mut condition_b =   
        &index_birth                           <   birth                        && 
        chx.key_2_filtration(&index_death)    ==  chx.key_2_filtration(&death);
        
        if condition_a || condition_b{
            // println!("index birth to A {:?} ",index_birth);
            // println!("index death to A {:?}",index_death);

            A.push(factored_complex.get_matched_basis_vector(1, &index_birth));
        }

    }

    let edge_size = good_edges.len();
    let triangle_size = good_triangles.len();
    let column_size_of_a = A.len(); 
    

    // Set program type
    let  mut program_type = Integer;
    if is_int{
        program_type = Integer;
    }
    else{
        program_type = Continuous;
    }

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

    // initialize the vector: x+
    let mut x_pos = Vec::new();
        
    for i in 0..edge_size{

        let  name = format!("{}{}", "x_pos", i);

        let  str_name = &name[..];
        // Should put weight in objective function because we don't want weight in constraint
        x_pos.push(model.add_var(str_name, program_type,1.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: x-
    let mut x_neg = Vec::new();
    
    for i in 0..edge_size{

        let  name = format!("{}{}", "x_neg", i);

        let  str_name = &name[..];
        x_neg.push(model.add_var(str_name, program_type, weight(index_2_edge.get(&i).unwrap()), 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: p 
    let mut p = Vec::new();
    
    for i in 0..column_size_of_a{

        let  name = format!("{}{}", "p", i);

        let  str_name = &name[..];
        p.push(model.add_var(str_name, program_type, 0.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // initialize the vector: q
    let mut q = Vec::new();
    
    for i in 0..triangle_size{

        let  name = format!("{}{}", "q", i);

        let  str_name = &name[..];
        q.push(model.add_var(str_name, program_type, 0.0, 0.0, INFINITY, &[], &[]).unwrap());
    }

    // Set objective function
    let mut obj_expression: LinExpr = LinExpr::new();
        
    for i in 0..edge_size{    
        
        obj_expression = obj_expression.add_term(1.0, x_pos[i].clone());

    }

    for i in 0..edge_size{    
        
        obj_expression = obj_expression.add_term(1.0, x_neg[i].clone());

    }
    
    // println!("obj expression: {:?}",obj_expression);


    // The cycle vector: x_orig
    let x_orig = factored_complex.get_matched_basis_vector(1, &birth);
    
   

    // Build Hashmap to record the index of the triangles
    let mut triangle_2_index: HashMap<MatrixIndexKey, usize> = HashMap::new();       
    let mut index_2_triangle: HashMap<usize, MatrixIndexKey> = HashMap::new();       
    // initialize indices to be 0   
    let mut tri_index:usize = 0;
    

    for triangle in good_triangles.iter() { // for each column
        if !triangle_2_index.contains_key(triangle) {
            triangle_2_index.insert(triangle.clone(), tri_index.clone());
            index_2_triangle.insert(tri_index.clone(), triangle.clone());
            tri_index = tri_index + 1;
        }
    }

    let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row, exhact::chx::ChxTransformKind::Boundary);

    // println!("good edges : {:?}",good_edges.clone());
    println!("good triangles : {:?}",good_triangles.clone());
    
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
        // println!("{:?} \n",constraint);
        model.add_constr("constraint", constraint, Equal, 0.0);
    }

    // Set up the objective function and let gurobi solve for the result.
    model.update().unwrap();
    model.set_objective(obj_expression,Minimize).unwrap();
    model.write("logfile.lp").unwrap();
    model.optimize().unwrap(); // Solve result
    let x_neg_val = model.get_values(attr::X, &x_neg).unwrap(); // Get the result for x+ and x-
    let x_pos_val = model.get_values(attr::X, &x_pos).unwrap();

    // println!("x_neg_val :{:?}",x_neg_val);
    // println!("x_pos_val :{:?}",x_pos_val);

    // Create a hashmap to store the answer and return 
    let mut ans: HashMap<MatrixIndexKey, f64> = HashMap::new();
    for i in 0..edge_size{
        if x_pos_val[i]- x_neg_val[i]!=0.{
            let edge_ans = index_2_edge.get(&i).unwrap().clone();
            let edge_val = x_pos_val[i]- x_neg_val[i];
            ans.insert(edge_ans, edge_val);
        }
    }

    return ans;

}



fn main() {  
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
    
    
    for j in 0..simplex_bar.len(){
        let birth = &simplex_bar[j].0;
        let death = &simplex_bar[j].1;
        // println!("birth: {:?} death: {:?}",birth,death);
        // Write solution to npy
    
        let solution_hash_edge = edge_opt(&factored_complex, birth,death, 1,false, |x| getLength(x, &dismat));
        let mut vertices_sol_vec = Vec::new();
        let mut coeff_sol_vec = Vec::new();
        // println!("weight {:?}",solution_hash_edge);

        for (print_key, print_val) in solution_hash_edge.iter() {
            vertices_sol_vec.push(print_key.vertices[0]);
            vertices_sol_vec.push(print_key.vertices[1]);
            coeff_sol_vec.push(*print_val);
        }
        let vertices_sol_arr = Array::from_vec(vertices_sol_vec);
        let coeff_sol_arr = Array::from_vec(coeff_sol_vec);
        let folder_name = format!("{}{}", "npy_files_dis_mat_cycle_", j);  
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

