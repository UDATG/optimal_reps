




pub enum SimplexWeights {
    Uniform,
    Volume,
}


pub fn tri_opt( stuff, weight: SimplexWeights )

{


    match weight {
        SimplexWeights::Uniform => {
            let weight_vector = ... make the weight vector
        }
        SimplexWeights::Volume => {
            let weight_vector = ... make the weight vector
        }
    }
}


x = tri_opt( stuff, SimplexWeights::Uniform)
