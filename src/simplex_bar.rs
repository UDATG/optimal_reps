use exhact::chx::{FactoredComplexBlockCsm, ChainComplex, Indexing};
use std::hash::Hash;
use std::fmt::Debug;
use std::ops::{Neg, AddAssign, Mul};


pub fn simplex_barcode< Filtration, MatrixIndexKey, SnzVal, OriginalChx >

    (   factored: & FactoredComplexBlockCsm< MatrixIndexKey, SnzVal, Filtration, OriginalChx >, 
        h_degree: usize  ) 
    -> 
    Vec<(MatrixIndexKey, MatrixIndexKey)> 
   
    where 
        Filtration: PartialOrd + Clone,
        MatrixIndexKey: PartialEq + Eq + Ord + Hash + Clone + Debug,
        SnzVal: Clone + PartialEq + Debug,
        OriginalChx: ChainComplex<MatrixIndexKey, SnzVal, Filtration>

{
    // let mut barcode = Vec::new();
    let mut simplexwise_barcode = Vec::new();

    if h_degree >= factored.dim_rowoper.len() { 
        return simplexwise_barcode; }
    //let keys = factored.original_complex.keys_ordered(h_degree);
    let indexing = &factored.dim_indexing[h_degree];
    let mut indexing_upper: &Indexing<MatrixIndexKey, MatrixIndexKey> = &Indexing::new();
    if h_degree+1 < factored.dim_rowoper.len() {
        indexing_upper = &factored.dim_indexing[h_degree+1];
    }

    


    for key in factored.original_complex.keys_unordered_itr(h_degree) { // for each k simplex
       
        if indexing.minkey_2_index.contains_key(&key) { 
            continue; 
        } // if this is a pivot column in d1, continue
        else if indexing_upper.majkey_2_index.contains_key(&key) {  // if this is a pivot row in d2
            let ind = indexing_upper.majkey_2_index[&key];          // d2 row index 
            let matched_key = &indexing_upper.index_2_minkey[ind];  // d2 column index 
            let diam1 = factored.original_complex.key_2_filtration(&key); // d1 birth time 
            let diam2 = factored.original_complex.key_2_filtration(&matched_key); // d2 death time 
            if diam1 == diam2 { 
                continue; 
            } // if equal, not a feature
            simplexwise_barcode.push((key.clone(), matched_key.clone()));   // else is a feature
         } 
         else {
            let diam1 = factored.original_complex.key_2_filtration(&key); // born 
            let diam2 = factored.original_complex.max_filtration(); // never dies 
            if diam1 == diam2 { 
                continue; 
            } 
            simplexwise_barcode.push((key.clone(), key.clone())); 
        }
        
    }
    return simplexwise_barcode;
}

