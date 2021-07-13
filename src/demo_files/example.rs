
fn main() {    

    // Step1 ; read distance matrix
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

    // Step2: build and factor the filtered chain complex
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
    
    

    // Step 3 optimize for a chosen feature
    for j in 2..3{//simplex_bar.len(){
        let birth = &simplex_bar[j].0;
        //println!("{:?}",factored_complex.get_matched_basis_vector(1, birth));
        

        let death = &simplex_bar[j].1;

        //uniform weight
        println!("uniform weight");
        let solution_hash_edge = tri_opt(true, 1, |x| 1., &factored_complex, birth, death);
        println!("Solution");
        for (print_key, print_val) in solution_hash_edge.iter() {
            println!("{:?}" ,(print_key, print_val));
        }
        
        // weight by area
        println!("weight by area");
        let solution_hash_edge = tri_opt(true, 1, |x| getArea(x, &dismat), &factored_complex, birth, death);
        println!("Solution");
        for (print_key, print_val) in solution_hash_edge.iter() {
            println!("{:?}" ,(print_key, print_val));
        }
    
    }


    // Alternative Step 3: optimize for a list of features
    let mut basis_unif = Vec::new();
    let mut basis_nonUnif = Vec::new();
    for j in 1..simplex_bar.len(){
        let birth = &simplex_bar[j].0;
       

        let death = &simplex_bar[j].1;

        //uniform weight
        println!("uniform weight");
        let solution_hash_edge = tri_opt(true, 1, |x| 1., &factored_complex, birth, death);
        basis.unif.push(solution_hash_edge);
        
        // weight by area
        println!("weight by area");
        let solution_hash_edge = tri_opt(true, 1, |x| getArea(x, &dismat), &factored_complex, birth, death);
        basis.nonUnif.push(solution_hash_edge);
    
    }

}