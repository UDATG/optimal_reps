// pub fn simplex_barcode(&self, h_degree: usize) -> Vec<(MatrixIndexKey, MatrixIndexKey)> {

//     // let mut barcode = Vec::new();
//     let mut simplexwise_barcode = Vec::new();

//     if h_degree >= self.dim_rowoper.len() { return simplexwise_barcode; }
//     //let keys = self.original_complex.keys_ordered(h_degree);
//     let indexing = &self.dim_indexing[h_degree];
//     let mut indexing_upper: &Indexing<MatrixIndexKey, MatrixIndexKey> = &Indexing::new();
//     if h_degree+1 < self.dim_rowoper.len() {
//         indexing_upper = &self.dim_indexing[h_degree+1];
//     }

//     for key in self.original_complex.keys_unordered_itr(h_degree) { // for each k simplex
//         if indexing.minkey_2_index.contains_key(&key) { continue; } // if this is a pivot column in d1, continue
//         else if indexing_upper.majkey_2_index.contains_key(&key) {  // if this is a pivot row in d2
//             let ind = indexing_upper.majkey_2_index[&key];          // d2 row index 
//             let matched_key = &indexing_upper.index_2_minkey[ind];  // d2 column index 
//             let diam1 = self.original_complex.key_2_filtration(&key); // d1 birth time 
//             let diam2 = self.original_complex.key_2_filtration(matched_key); // d2 death time 
//             if diam1 == diam2 { continue; } // if equal, not a feature
//             simplexwise_barcode.push((key.clone(), matched_key.clone()));   // else is a feature
//          } 
//          else {
//             let diam1 = self.original_complex.key_2_filtration(&key); // born 
//             let diam2 = self.original_complex.max_filtration(); // never dies 
//             if diam1 == diam2 { continue; } 
//             simplexwise_barcode.push((key.clone(), key.clone())); 
//         }
//     }
//     return simplexwise_barcode;
// }